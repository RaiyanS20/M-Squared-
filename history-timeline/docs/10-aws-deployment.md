# 10 — Deploying to AWS

> **Concept** → The code → Do it yourself → Check yourself

This chapter is a guide, not a script that has been run against your account. The
artifacts in `infra/` and `server/Dockerfile` are real and complete; the AWS
resources are yours to create. **Set a budget alarm before you start** — the
section on cost is not an afterthought.

## The target

```
                    ┌──────────────┐
   Browser ────────▶│  CloudFront  │  CDN, TLS, one domain
                    └──┬────────┬──┘
             /* (static)│        │ /api/* (dynamic)
                  ┌─────▼───┐  ┌─▼──────────────────┐
                  │   S3    │  │ Application Load    │
                  │ (React) │  │ Balancer            │
                  └─────────┘  └─┬──────────────────┘
                                 │
                          ┌──────▼──────────┐
                          │ ECS Fargate     │  Node.js API, 2+ tasks
                          └──────┬──────────┘
                                 │
                          ┌──────▼──────────┐
                          │ RDS PostgreSQL  │  private subnet
                          └─────────────────┘
```

**Why one CloudFront distribution for both?** The browser sees a single domain,
so there is no CORS in production and no preflight request on every call. It also
means one TLS certificate and one cache to reason about.

## Why these services

| Need | Service | Why not the alternative |
|---|---|---|
| Static React files | **S3 + CloudFront** | A server to send unchanging files is waste. The CDN puts them near the user. |
| Run the API | **ECS Fargate** | No servers to patch. Lambda is cheaper at low traffic but adds cold starts and a different programming model — and here you would be re-learning the deployment, not the app. |
| PostgreSQL | **RDS** | Managed backups, patching and failover. Running Postgres on EC2 means you are now a DBA. |
| Secrets | **Secrets Manager** | Rotatable, audited, never in an image or a task definition. |
| Logs | **CloudWatch Logs** | ECS ships them with two lines of config. |

**Fargate vs Lambda** is the interesting trade-off. Lambda would be nearly free
at this traffic, but our Express app assumes a long-lived process and a
connection pool — neither fits Lambda's model well. Fargate runs the code you
already have. If this were a fresh greenfield API with spiky traffic, Lambda
would be the better answer.

## Step 1 — Build and push the API image

`server/Dockerfile` is a multi-stage build: stage 1 compiles TypeScript, stage 2
takes only `dist/` plus production dependencies. No compiler in production, a
much smaller image, `USER node` so nothing runs as root.

```bash
aws ecr create-repository --repository-name history-timeline-api

aws ecr get-login-password --region "$REGION" \
  | docker login --username AWS --password-stdin "$ACCOUNT.dkr.ecr.$REGION.amazonaws.com"

docker build -f server/Dockerfile -t history-timeline-api:v1 .
docker tag history-timeline-api:v1 "$ACCOUNT.dkr.ecr.$REGION.amazonaws.com/history-timeline-api:v1"
docker push "$ACCOUNT.dkr.ecr.$REGION.amazonaws.com/history-timeline-api:v1"
```

**Tag with a version or a git SHA, never `latest`.** `latest` means you cannot
say what is running, and cannot roll back to what was.

## Step 2 — RDS PostgreSQL

```bash
aws rds create-db-instance \
  --db-instance-identifier history-timeline \
  --db-instance-class db.t4g.micro \
  --engine postgres --engine-version 16 \
  --allocated-storage 20 \
  --master-username historian \
  --manage-master-user-password \
  --no-publicly-accessible \
  --backup-retention-period 7
```

Non-negotiables:

- **`--no-publicly-accessible`.** A database on the public internet is scanned
  within minutes. There is no configuration that makes this acceptable.
- **Security group allows 5432 only from the ECS task's security group** —
  reference the group, not a CIDR.
- **`--manage-master-user-password`** puts the password in Secrets Manager
  instead of your shell history.
- **Backups on**, with a retention period you chose deliberately.

Then store the connection string:

```bash
aws secretsmanager create-secret --name history-timeline/database-url \
  --secret-string "postgres://historian:PASSWORD@ENDPOINT:5432/history_timeline"
```

## Step 3 — Run migrations

Migrations run **before** the new code, as a one-off ECS task using the same
image:

```bash
aws ecs run-task --cluster history-timeline --launch-type FARGATE \
  --task-definition history-timeline-api \
  --overrides '{"containerOverrides":[{"name":"api","command":["node","server/dist/scripts/migrate.js"]}]}' \
  --network-configuration "awsvpcConfiguration={subnets=[$SUBNETS],securityGroups=[$SG]}"
```

Same image, different command — so the migration runs against exactly the code
being deployed.

**The rule that makes rolling deploys safe:** during a rollout, old and new code
run *simultaneously*. Every migration must be backwards compatible with the
currently running version.

- ✅ Add a nullable column
- ✅ Add a table, add an index
- ❌ Drop or rename a column in the same deploy that stops using it

Renaming is three deploys: add the new column and write to both → backfill and
read from the new one → drop the old. Tedious, and the alternative is downtime.

## Step 4 — ECS Fargate

`infra/task-definition.json` is ready to fill in. The two details that matter:

```jsonc
"secrets": [{ "name": "DATABASE_URL", "valueFrom": "arn:aws:secretsmanager:..." }]
```

**`secrets`, not `environment`.** Values in `environment` are visible to anyone
who can describe the task. `secrets` are fetched at start-up and never appear in
the definition.

```jsonc
"cpu": "256", "memory": "512"
```

Per task, and Fargate only accepts specific pairings. This is the smallest, and
is ample here.

Create the service with at least two tasks across two availability zones:

```bash
aws ecs create-service --cluster history-timeline \
  --service-name api --task-definition history-timeline-api \
  --desired-count 2 --launch-type FARGATE \
  --health-check-grace-period-seconds 30 \
  --load-balancers "targetGroupArn=$TG_ARN,containerName=api,containerPort=4000"
```

**Point the ALB health check at `/ready`, not `/health`.** `/ready` checks the
database, so a task that cannot reach RDS is pulled out of rotation instead of
serving 500s. That distinction (chapter 6) exists for exactly this line of
configuration.

### Connection pool maths

`PostgresDatabase` uses `max: 10` **per process**. Two tasks = 20 connections. A
`db.t4g.micro` allows roughly 80. Scale to twenty tasks without thinking and you
exhaust connections before CPU — an outage that looks like nothing is wrong.

Either lower `max`, size the instance for `tasks × max`, or put RDS Proxy in
front. **Know which one you are relying on.**

### Graceful shutdown

`main.ts` handles SIGTERM: stop accepting connections, finish in-flight requests,
close the pool, exit.

```ts
process.on('SIGTERM', () => void shutdown('SIGTERM'));
```

Without it, ECS kills the process mid-request on every deploy and users see
connection resets. `stopTimeout: 30` in the task definition gives it room.

## Step 5 — The React app

```bash
export S3_BUCKET=my-history-timeline-site
export CLOUDFRONT_DISTRIBUTION_ID=E1234567890ABC
./infra/deploy-client.sh
```

The script uploads in two passes, and the reason is worth internalising:

```bash
# Hashed filenames (index-a1b2c3.js) — a change produces a NEW name
--cache-control "public, max-age=31536000, immutable"

# index.html — points at the hashes, so it must never be cached
--cache-control "no-cache, no-store, must-revalidate"
```

Cache `index.html` and users keep loading the old app while the new assets sit
unused. This is the single most common static-deploy bug.

### CloudFront behaviours

| Path pattern | Origin | Cache |
|---|---|---|
| `/api/*` | ALB | disabled — forward all headers |
| `/*` | S3 | default |

Also set a **custom error response**: 404 → `/index.html` with status 200, so
client-side routes work if you add a router later.

## Step 6 — Watch it

The four alarms worth having on day one:

| Alarm | Threshold | Means |
|---|---|---|
| ALB 5xx rate | > 1% for 5 min | The API is failing |
| ECS CPU | > 80% for 10 min | Scale out |
| RDS free storage | < 10% | Act before it stops accepting writes |
| RDS connections | > 80% of max | The pool maths above is wrong |

Query logs with CloudWatch Logs Insights:

```
fields @timestamp, @message
| filter @message like /unhandled error/
| sort @timestamp desc | limit 50
```

## What this costs

Rough monthly figures for a low-traffic educational site (us-east-1, 2025 list
prices — check current pricing):

| Service | Configuration | ~USD/month |
|---|---|---|
| ECS Fargate | 2 × 0.25 vCPU / 0.5 GB | ~15 |
| RDS PostgreSQL | db.t4g.micro, 20 GB, single-AZ | ~15 |
| ALB | one, low traffic | ~18 |
| S3 + CloudFront | a few GB out | ~2 |
| Secrets Manager | one secret | ~0.40 |
| **Total** | | **~50** |

The **ALB is the surprise** — a fixed ~$18/month whether or not anyone visits.
For a learning project, alternatives: run one Fargate task with a public IP and
no ALB (no zero-downtime deploys), or use App Runner, or put the API on Lambda +
API Gateway and pay per request.

Multi-AZ RDS roughly doubles the database cost and is right for production, hard
to justify for a course project.

**Do this now, before anything else:**

```bash
aws budgets create-budget --account-id "$ACCOUNT" --budget \
  '{"BudgetName":"history-timeline","BudgetLimit":{"Amount":"25","Unit":"USD"},
    "TimeUnit":"MONTHLY","BudgetType":"COST"}'
```

Then delete what you are not using. `aws ecs update-service --desired-count 0`
and stopping the RDS instance (it restarts automatically after 7 days) will cut
most of it. **The most common AWS learning experience is a surprise bill.**

## Security checklist

- [ ] RDS not publicly accessible; SG references the ECS SG, not a CIDR
- [ ] Secrets in Secrets Manager, never in the image or task definition
- [ ] Container runs as `USER node`
- [ ] `CORS_ORIGINS` names your domain — never `*`
- [ ] TLS everywhere: ACM certificate on CloudFront, `DATABASE_SSL=true`
- [ ] IAM roles scoped to what the task actually needs
- [ ] Stack traces never returned in production (`isProduction` in `errorHandler`)
- [ ] Budget alarm set

## Do it yourself

1. **Build and run the image locally** — no AWS account needed:
   ```bash
   docker build -f server/Dockerfile -t history-timeline-api .
   docker run --rm -p 4000:4000 \
     -e DATABASE_URL='postgres://historian:historian@host.docker.internal:5432/history_timeline' \
     history-timeline-api
   curl localhost:4000/ready
   ```
   Compare `docker images` sizes for the build and runtime stages.

2. **Test graceful shutdown.** Start the container, `docker stop` it, and watch
   the log line. Then remove the SIGTERM handler and compare.

3. **Deploy the client only.** S3 + CloudFront is the cheapest useful half. Point
   `VITE_API_BASE_URL` at an API running anywhere.

4. **Plan a column rename** across three deploys, writing out what runs at each
   step. This is the exercise that makes zero-downtime deployment click.

## Check yourself

- Why must the ALB health check use `/ready` rather than `/health`?
- Why must migrations be backwards compatible?
- Why is `latest` a bad image tag?

<details>
<summary>Answers</summary>

- **`/ready`**: it checks the database. `/health` only says the process is
  alive, so a task that cannot reach RDS would stay in rotation returning 500s.
  `/ready` returns 503 and the ALB stops sending it traffic.
- **Backwards compatible migrations**: during a rolling deploy both versions run
  at once against one database. A migration the old code cannot survive breaks
  every request still being served by an old task.
- **`latest`**: it is mutable. You cannot say what is running, you cannot roll
  back to a known-good build, and two tasks started minutes apart can be running
  different code with the same tag.

</details>

---

## You have reached the end

What you built: a full-stack, tested, deployable web application with a layered
architecture, a real domain model, 132 tests and a deployment path.

What to do next, in rough order of value:

1. **Widen the timeline to 100 AD.** Change `Year.MVP_START`, add pre-1900
   events, and solve the UI problem 1,900 year-buttons creates (zoom levels? a
   century rail that expands into decades?). This is a genuinely interesting
   design problem and the architecture is ready for it.
2. **Add a router** and deep links — `/1994/ZAF` should be shareable.
3. **Add search** across event titles and summaries (PostgreSQL full-text search
   is right here; `tsvector` and a GIN index).
4. **Let people contribute events**, with the source URL required — which means
   authentication, moderation, and your first `POST` endpoint.
5. **Deploy it** and put the URL somewhere.

Each of those touches a different layer, and you now know which one.

← Back to [Start here](00-start-here.md)
