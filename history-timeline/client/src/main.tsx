import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { App } from './App.js';

/**
 * The client's entry point — the frontend equivalent of the server's
 * composition root: find the mount node, render the tree, and nothing else.
 *
 * StrictMode intentionally double-invokes effects in development to surface
 * missing cleanup. Our `useAsyncResource` aborts its request on cleanup, so it
 * passes that check — which is exactly the point of the exercise.
 */
const container = document.getElementById('root');
if (!container) {
  throw new Error('Root element #root is missing from index.html');
}

createRoot(container).render(
  <StrictMode>
    <App />
  </StrictMode>,
);
