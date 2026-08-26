// `tsc` compiles .ts and ignores everything else, so the .sql migrations would
// never reach dist/. This copies them after a build.
import { cp } from 'node:fs/promises';

await cp('src/db/migrations', 'dist/db/migrations', { recursive: true });
console.log('[build] copied migrations to dist/db/migrations');
