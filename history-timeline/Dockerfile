# One image that serves both the API and the browser app.
#
# There is no build stage, because there is nothing to build — no bundler, no
# transpiler. That is a genuine benefit of the vanilla stack: the code that runs
# in production is the code you wrote and read.

FROM node:22-alpine
WORKDIR /app
ENV NODE_ENV=production

# Copy manifests first. Docker caches layers, so dependencies are only
# reinstalled when the manifests change — not on every source edit.
COPY package.json package-lock.json ./
RUN npm ci --omit=dev && npm cache clean --force

COPY server ./server
COPY client ./client

# Never run as root. If the process is compromised, this is the difference
# between "an attacker has the app" and "an attacker has the container".
USER node

EXPOSE 4000

# Fails the container health check if the API stops answering, so the platform
# restarts it.
HEALTHCHECK --interval=30s --timeout=3s --start-period=10s --retries=3 \
  CMD node -e "fetch('http://localhost:4000/health').then(r=>process.exit(r.ok?0:1)).catch(()=>process.exit(1))"

CMD ["node", "server/src/main.js"]
