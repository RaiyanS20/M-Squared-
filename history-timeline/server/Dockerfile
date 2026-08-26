# Multi-stage build for the API.
#
# Stage 1 has the full toolchain (TypeScript, dev dependencies). Stage 2 gets
# only the compiled output and production dependencies — a smaller image, a
# smaller attack surface, and no compiler shipped to production.

# ---- Stage 1: build -------------------------------------------------------
FROM node:22-alpine AS build
WORKDIR /app

# Copy manifests first. Docker caches layers, so dependencies are only
# reinstalled when the manifests change — not on every source edit.
COPY package.json package-lock.json ./
COPY server/package.json ./server/
COPY client/package.json ./client/
RUN npm ci --workspace=server --include-workspace-root

COPY server ./server
RUN npm run build --workspace=server

# ---- Stage 2: runtime -----------------------------------------------------
FROM node:22-alpine AS runtime
WORKDIR /app
ENV NODE_ENV=production

# Production dependencies only: no TypeScript, no vitest, no supertest.
COPY package.json package-lock.json ./
COPY server/package.json ./server/
RUN npm ci --omit=dev --workspace=server --include-workspace-root && npm cache clean --force

COPY --from=build /app/server/dist ./server/dist

# Never run as root. If the process is compromised, this is the difference
# between "an attacker has the app" and "an attacker has the container".
USER node

EXPOSE 4000

# Fails the container health check if the API stops answering. ECS restarts it.
HEALTHCHECK --interval=30s --timeout=3s --start-period=10s --retries=3 \
  CMD node -e "fetch('http://localhost:4000/health').then(r=>process.exit(r.ok?0:1)).catch(()=>process.exit(1))"

CMD ["node", "server/dist/main.js"]
