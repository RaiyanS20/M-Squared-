#!/usr/bin/env bash
#
# Build the React app and publish it to S3 behind CloudFront.
#
# The two-pass upload is the important part. Hashed asset filenames
# (index-a1b2c3.js) can be cached for a year because a change produces a NEW
# filename. index.html must never be cached, because it is the file that points
# at the new hashes — cache it and users keep loading the old app.
set -euo pipefail

BUCKET="${S3_BUCKET:?set S3_BUCKET, e.g. my-history-timeline-site}"
DISTRIBUTION="${CLOUDFRONT_DISTRIBUTION_ID:?set CLOUDFRONT_DISTRIBUTION_ID}"

cd "$(dirname "$0")/.."

echo "==> Building the client"
npm run build --workspace=client

echo "==> Uploading hashed assets (immutable, 1 year)"
aws s3 sync client/dist "s3://${BUCKET}" \
  --delete \
  --exclude "index.html" \
  --cache-control "public, max-age=31536000, immutable"

echo "==> Uploading index.html (never cached)"
aws s3 cp client/dist/index.html "s3://${BUCKET}/index.html" \
  --cache-control "no-cache, no-store, must-revalidate"

echo "==> Invalidating CloudFront"
aws cloudfront create-invalidation \
  --distribution-id "${DISTRIBUTION}" \
  --paths "/index.html" \
  --query 'Invalidation.Id' --output text

echo "==> Done"
