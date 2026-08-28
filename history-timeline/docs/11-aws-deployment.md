# 11 — Deploying to AWS

> **Concept** → The code → Do it yourself → Check yourself

This chapter is a guide, not a script that has been run against your account. The
artifacts in `infra/` and the `Dockerfile` are real and complete; the AWS
resources are yours to create. **Set a budget alarm before you start** — the
section on cost is not an afterthought.

## Two shapes, and start with the simple one

**Simple — one container serves everything.** Exactly as in development: the Node
process serves the API and the static files.

```
Browser → ALB → ECS Fargate (API + static files) → RDS PostgreSQL
```

One thing to deploy. No CORS. No CDN to invalidate. This is a completely
respectable production setup for a site like this.

**Split — static files on a CDN.** `client/` goes to S3 behind CloudFront; the
API stays on Fargate. Faster worldwide, cheaper at scale, and the arrangement you
will meet in industry. It is also where `CORS_ORIGINS` finally matters, because
the page and the API are now on different origins.

**Do the simple one first.** Move to the split one when you can say why you are
moving — that sentence is the actual learning.

## Why these services

| Need | Service | Why not the alternative |
|---|---|---|
| Run the app | **ECS Fargate** | No servers to patch. Lambda is cheaper at low traffic but adds cold starts and does not fit a long-lived connection pool. |
| PostgreSQL | **RDS** | Managed backups, patching and failover. Running Postgres on EC2 means you are now a DBA. |
| Secrets | **Secrets Manager** | Rotatable, audited, never in an image. |
| Logs | **CloudWatch Logs** | ECS ships them with two lines of config. |
| Static files (split) | **S3 + CloudFront** | A server to send unchanging files is waste. |

**Fargate vs Lambda** is the interesting trade-off. Lambda would be nearly free
at this traffic, but our Express app assumes a long-lived process and a
connection pool — neither fits Lambda's model well. Fargate runs the code you
already have.

## Step 1 — Build and push the image

There is **no build stage** in the Dockerfile, because there is nothing to
build. The code that runs in production is the code you wrote and read — a real
benefit of the vanilla stack.

```bash
aws ecr create-repository --repository-name history-timeline

aws ecr get-login-password --region "$REGION" \
  | docker login --username AWS --password-stdin "$ACCOUNT.dkr.ecr.$REGION.amazonaws.com"

docker build -t history-timeline:v1 .
docker tag history-timeline:v1 "$ACCOUNT.dkr.ecr.$REGION.amazonaws.com/history-timeline:v1"
docker push "$ACCOUNT.dkr.ecr.$REGION.amazonaws.com/history-timeline:v1"
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

## Step 3 — Migrate and seed

Migrations run **before** the new code, as a one-off ECS task using the same
image:

```bash
aws ecs run-task --cluster history-timeline --launch-type FARGATE \
  --task-definition history-timeline \
  --overrides '{"containerOverrides":[{"name":"app","command":["node","server/src/scripts/migrate.js"]}]}' \
  --network-configuration "awsvpcConfiguration={subnets=[$SUBNETS],securityGroups=[$SG]}"
```

Same image, different command — so the migration runs against exactly the code
being deployed. Run `seed.js` the same way; it is idempotent, so re-running it is
safe.

**The rule that makes rolling deploys safe:** during a rollout, old and new code
run *simultaneously*. Every migration must be backwards compatible with the
version currently running.

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
who can describe the task.

```jsonc
"cpu": "256", "memory": "512"
```

Per task, and Fargate accepts only specific pairings. This is the smallest, and
is ample here.

Create the service with at least two tasks across two availability zones:

```bash
aws ecs create-service --cluster history-timeline \
  --service-name app --task-definition history-timeline \
  --desired-count 2 --launch-type FARGATE \
  --health-check-grace-period-seconds 30 \
  --load-balancers "targetGroupArn=$TG_ARN,containerName=app,containerPort=4000"
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

`main.js` handles SIGTERM: stop accepting connections, finish in-flight requests,
close the pool, exit. Without it, ECS kills the process mid-request on every
deploy and users see connection resets. `stopTimeout: 30` gives it room.

## Step 5 — (Split shape only) the browser app on a CDN

```bash
export S3_BUCKET=my-history-timeline-site
export CLOUDFRONT_DISTRIBUTION_ID=E1234567890ABC
./infra/deploy-client.sh
```

The script uploads in two passes, and the reason is worth internalising:

```bash
# css/ and js/ — short cache, revalidated
--cache-control "public, max-age=300, must-revalidate"

# index.html — never cached; it points at everything else
--cache-control "no-cache, no-store, must-revalidate"
```

Cache `index.html` aggressively and users keep loading the old page. This is the
single most common static-deploy bug.

> **A note on cache busting.** Because this project does not hash filenames, the
> assets get a short cache rather than a year-long one. If you later add a build
> step that emits `app-a1b2c3.js`, you can cache those for a year with
> `immutable`, because a change produces a new *name*. That is the main thing a
> bundler buys you here, and it is worth knowing what you are trading.

Then set two CloudFront behaviours and switch the API's `CORS_ORIGINS` to your
domain:

| Path pattern | Origin | Cache |
|---|---|---|
| `/api/*` | ALB | disabled — forward all headers |
| `/*` | S3 | default |

## Step 6 — Watch it

The four alarms worth having on day one:

| Alarm | Threshold | Means |
|---|---|---|
| ALB 5xx rate | > 1% for 5 min | The app is failing |
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
| Secrets Manager | one secret | ~0.40 |
| **Simple shape total** | | **~48** |
| S3 + CloudFront | a few GB out (split shape) | ~2 |

The **ALB is the surprise** — a fixed ~$18/month whether or not anyone visits.
For a learning project the alternatives are: one Fargate task with a public IP
and no ALB (no zero-downtime deploys), or AWS App Runner, or a small VPS
elsewhere for a few dollars.

Multi-AZ RDS roughly doubles the database cost and is right for production, hard
to justify for a course project.

**Do this now, before anything else:**

```bash
aws budgets create-budget --account-id "$ACCOUNT" --budget \
  '{"BudgetName":"history-timeline","BudgetLimit":{"Amount":"25","Unit":"USD"},
    "TimeUnit":"MONTHLY","BudgetType":"COST"}'
```

Then delete what you are not using. `aws ecs update-service --desired-count 0`
and stopping the RDS instance (it restarts automatically after 7 days) cuts most
of it. **The most common AWS learning experience is a surprise bill.**

## Security checklist

- [ ] RDS not publicly accessible; SG references the ECS SG, not a CIDR
- [ ] Secrets in Secrets Manager, never in the image or task definition
- [ ] Container runs as `USER node`
- [ ] `CORS_ORIGINS` names your domain — never `*` (split shape only)
- [ ] TLS everywhere: ACM certificate on the ALB/CloudFront, `DATABASE_SSL=true`
- [ ] IAM roles scoped to what the task actually needs
- [ ] Stack traces never returned in production (`isProduction` in `errorHandler`)
- [ ] Budget alarm set

## Do it yourself

1. **Build and run the image locally** — no AWS account needed:
   ```bash
   docker build -t history-timeline .
   docker run --rm -p 4000:4000 \
     -e DATABASE_URL='postgres://historian:historian@host.docker.internal:5432/history_timeline' \
     history-timeline
   curl localhost:4000/ready
   ```
   Then open <http://localhost:4000> — the container serves the page too.

2. **Test graceful shutdown.** Start the container, `docker stop` it, and watch
   the log line. Then remove the SIGTERM handler and compare.

3. **Plan a column rename** across three deploys, writing out what runs at each
   step. This is the exercise that makes zero-downtime deployment click.

4. **Price it yourself** in the AWS pricing calculator for your region, and
   decide whether you would run the simple or split shape.

## Check yourself

- Why must the ALB health check use `/ready` rather than `/health`?
- Why must migrations be backwards compatible?
- Why is `latest` a bad image tag?

<details>
<summary>Answers</summary>

- **`/ready`**: it checks the database. `/health` only says the process is alive,
  so a task that cannot reach RDS would stay in rotation returning 500s.
- **Backwards compatible migrations**: during a rolling deploy both versions run
  at once against one database. A migration the old code cannot survive breaks
  every request still being served by an old task.
- **`latest`**: it is mutable. You cannot say what is running, cannot roll back to
  a known-good build, and two tasks started minutes apart can be running
  different code under the same tag.

</details>

---

## You have reached the end

What you built: a full-stack, tested, deployable web application — vanilla
JavaScript, HTML and CSS in the browser, a layered OOP Node.js API, PostgreSQL,
163 tests, and a deployment path. No build step, five dependencies.

What to do next, in rough order of value:

1. **Widen the timeline to 100 AD.** Change `Year.MVP_START`, add pre-1900
   events, and solve the UI problem 1,900 year-buttons creates (zoom levels? a
   century rail that expands into decades?). A genuinely interesting design
   problem, and the architecture is ready for it.
2. **Add deep links** — `/1994/ZAF` should be shareable, using `history.pushState`
   and `popstate`. This is where a router starts to earn its place.
3. **Add search** across titles and summaries (PostgreSQL full-text search with
   `tsvector` and a GIN index is exactly right here).
4. **Let people contribute events**, with the source URL required — which means
   authentication, moderation, and your first `POST` endpoint.
5. **Add ESLint, then TypeScript.** In that order. You will appreciate what
   types give you far more having felt their absence, and the architecture is
   already shaped for them.
6. **Deploy it** and put the URL somewhere.

Each of those touches a different layer, and you now know which one.

← Back to [Start here](00-start-here.md)
