#!/usr/bin/env bash
#
# Publish the browser app to S3 behind CloudFront (the "split" shape).
#
# There is no build step — the files in client/ are the files that ship. That is
# a genuine benefit of the vanilla stack.
#
# The two-pass upload is the important part. Hashed or versioned assets could be
# cached for a year; index.html must NEVER be cached, because it is the file
# that points at everything else. Since this project does not hash filenames,
# both CSS and JS are given a short cache with revalidation.
set -euo pipefail

BUCKET="${S3_BUCKET:?set S3_BUCKET, e.g. my-history-timeline-site}"
DISTRIBUTION="${CLOUDFRONT_DISTRIBUTION_ID:?set CLOUDFRONT_DISTRIBUTION_ID}"

cd "$(dirname "$0")/.."

echo "==> Uploading css/ and js/"
aws s3 sync client "s3://${BUCKET}" \
  --delete \
  --exclude "index.html" \
  --exclude "tests/*" \
  --cache-control "public, max-age=300, must-revalidate"

echo "==> Uploading index.html (never cached)"
aws s3 cp client/index.html "s3://${BUCKET}/index.html" \
  --cache-control "no-cache, no-store, must-revalidate"

echo "==> Invalidating CloudFront"
aws cloudfront create-invalidation \
  --distribution-id "${DISTRIBUTION}" \
  --paths "/index.html" "/js/*" "/css/*" \
  --query 'Invalidation.Id' --output text

echo "==> Done"
