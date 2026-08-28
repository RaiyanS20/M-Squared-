# Infrastructure

Deployment artifacts and notes for running this project on AWS. The narrative
walkthrough is [docs/11 — Deploying to AWS](../docs/11-aws-deployment.md); this
directory holds the concrete pieces.

## Contents

| File | What it is |
|---|---|
| `../Dockerfile` | One image serving both the API and the browser app |
| `task-definition.json` | ECS Fargate task definition (fill in the placeholders) |
| `deploy-client.sh` | Publish `client/` to S3 + CloudFront (the split option) |

## Two deployment shapes

**Simple — one container serves everything.** The Node process serves the API
*and* the static files, exactly as it does in development. One thing to deploy,
no CORS, no CDN to invalidate. Good enough for a learning project and a real
amount of traffic.

```
Browser → ALB → ECS Fargate (API + static) → RDS
```

**Split — static files on a CDN.** The browser app goes to S3 behind CloudFront;
the API stays on Fargate. Faster worldwide, cheaper at scale, and the setup you
will meet in industry.

```
                    ┌──────────────┐
   Browser ────────▶│  CloudFront  │
                    └──┬────────┬──┘
             /* (static)│        │ /api/*
                  ┌─────▼───┐  ┌─▼──────────────┐
                  │   S3    │  │      ALB       │
                  └─────────┘  └─┬──────────────┘
                                 │
                          ┌──────▼──────────┐
                          │ ECS Fargate     │
                          └──────┬──────────┘
                                 │
                          ┌──────▼──────────┐
                          │ RDS PostgreSQL  │  private subnet
                          └─────────────────┘
```

Start with the simple shape. Move to the split one when you can say why.

## Before you deploy

- [ ] `npm test` passes, including the Postgres contract tests
- [ ] `DATABASE_URL` is in Secrets Manager, **not** in the task definition
- [ ] The RDS instance is in a private subnet with no public IP
- [ ] The RDS security group allows 5432 **only** from the ECS task's security
      group — not from `0.0.0.0/0`
- [ ] The ALB health check points at `/ready`, not `/health`
- [ ] `CORS_ORIGINS` is set only if you chose the split shape
- [ ] Automated RDS backups are on, with a retention period you have chosen
- [ ] A budget alarm exists (see chapter 11 — this is the one people skip)

## Deploy order

Sequence matters. Migrations run **before** the new code:

1. `npm run db:migrate` against the production database (as a one-off ECS task)
2. Push the new image to ECR
3. Update the ECS service — it rolls tasks over gradually
4. (Split shape only) sync `client/` to S3 and invalidate CloudFront

This ordering only works if each migration is **backwards compatible** with the
currently running code — during a rolling deploy, old and new run side by side.
Adding a nullable column is safe. Dropping one is not: split it across two
deploys (stop using it, then drop it).
