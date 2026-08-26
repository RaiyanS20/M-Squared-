# Infrastructure

Deployment artifacts and notes for running this project on AWS. The narrative
walkthrough is [docs/10 — Deploying to AWS](../docs/10-aws-deployment.md); this
directory holds the concrete pieces.

## Contents

| File | What it is |
|---|---|
| `../server/Dockerfile` | Multi-stage build for the API image |
| `task-definition.json` | ECS Fargate task definition (fill in the placeholders) |
| `deploy-client.sh` | Build and publish the React app to S3 + CloudFront |

## The target architecture

```
                    ┌──────────────┐
   Browser ────────▶│  CloudFront  │  CDN, TLS, caching
                    └──┬────────┬──┘
             /* (static)│        │ /api/* (dynamic)
                  ┌─────▼───┐  ┌─▼──────────────────┐
                  │   S3    │  │ Application Load    │
                  │ (React) │  │ Balancer            │
                  └─────────┘  └─┬──────────────────┘
                                 │
                          ┌──────▼──────────┐
                          │ ECS Fargate     │  the Node.js API
                          │ (2+ tasks)      │
                          └──────┬──────────┘
                                 │
                          ┌──────▼──────────┐
                          │ RDS PostgreSQL  │  private subnet, no public IP
                          └─────────────────┘
```

One CloudFront distribution serves both origins, so the browser sees a single
domain and there is no CORS in production.

## Before you deploy

- [ ] `npm test` passes, including the Postgres contract tests
- [ ] `npm run build` succeeds in both workspaces
- [ ] `DATABASE_URL` is in Secrets Manager, **not** in the task definition
- [ ] The RDS instance is in a private subnet with no public IP
- [ ] The RDS security group allows 5432 **only** from the ECS task's security
      group — not from `0.0.0.0/0`
- [ ] The ALB health check points at `/ready`, not `/health`
- [ ] Automated RDS backups are on, with a retention period you have chosen
- [ ] A budget alarm exists (see chapter 10 — this is the one people skip)

## Deploy order

Sequence matters. Migrations run **before** the new code:

1. `npm run db:migrate` against the production database (as a one-off ECS task)
2. Push the new image to ECR
3. Update the ECS service — it rolls tasks over gradually
4. Sync the client build to S3 and invalidate the CloudFront cache

This ordering only works if each migration is **backwards compatible** with the
currently running code — during a rolling deploy, old and new run side by side.
Adding a nullable column is safe. Dropping a column is not: split it into two
deploys (stop using it, then drop it).
