import { el } from '../dom.js';

/**
 * Every remote request has the same three unhappy paths, so they get one module
 * instead of being re-invented in each panel.
 *
 * `role="status"` with `aria-live="polite"` makes a screen reader ANNOUNCE the
 * text when it appears. Without it, a blind user gets silence where a sighted
 * user sees a spinner.
 */

export function loadingMessage(label) {
  return el('p', { class: 'status status--loading', role: 'status', 'aria-live': 'polite' }, [
    el('span', { class: 'status__spinner', 'aria-hidden': 'true' }),
    label,
  ]);
}

export function errorMessage(error, onRetry) {
  return el('div', { class: 'status status--error', role: 'alert' }, [
    el('p', { class: 'status__title' }, error.message),
    // Retrying a 400 will fail identically, so only offer it when it might help.
    onRetry && !error.isClientError
      ? el('button', { type: 'button', class: 'button', on: { click: onRetry } }, 'Try again')
      : null,
  ]);
}

export function emptyMessage(text) {
  return el('p', { class: 'status status--empty', role: 'status' }, text);
}
