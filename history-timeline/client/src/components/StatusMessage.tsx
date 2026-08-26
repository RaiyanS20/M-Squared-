import type { ApiError } from '../api/index.js';

/**
 * Every remote request has the same three unhappy paths, so they get one
 * component instead of being re-invented in each screen.
 *
 * `role="status"` and `aria-live` mean a screen reader announces "Loading…" and
 * error text when it appears, rather than leaving the user staring at silence.
 */

export function Loading({ label }: { label: string }) {
  return (
    <p className="status status--loading" role="status" aria-live="polite">
      <span className="status__spinner" aria-hidden="true" />
      {label}
    </p>
  );
}

export function ErrorMessage({ error, onRetry }: { error: ApiError; onRetry?: () => void }) {
  return (
    <div className="status status--error" role="alert">
      <p className="status__title">{error.message}</p>
      {/* Retrying a 400 will fail identically, so only offer it when it might help. */}
      {onRetry && !error.isClientError && (
        <button type="button" className="button" onClick={onRetry}>
          Try again
        </button>
      )}
    </div>
  );
}

export function EmptyMessage({ children }: { children: React.ReactNode }) {
  return (
    <p className="status status--empty" role="status">
      {children}
    </p>
  );
}
