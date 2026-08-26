import { useEffect, useReducer } from 'react';
import { ApiError } from '../api/index.js';

/**
 * The three states any remote data can be in. Modelling them as a discriminated
 * union (rather than three loose booleans) makes the impossible states —
 * "loading AND has an error" — unrepresentable, and lets TypeScript narrow
 * `data` to non-null inside the success branch.
 */
export type AsyncState<T> =
  | { status: 'idle'; data: null; error: null }
  | { status: 'loading'; data: null; error: null }
  | { status: 'success'; data: T; error: null }
  | { status: 'error'; data: null; error: ApiError };

type Action<T> =
  | { type: 'idle' }
  | { type: 'loading' }
  | { type: 'success'; data: T }
  | { type: 'error'; error: ApiError };

const IDLE = { status: 'idle', data: null, error: null } as const;

function reducer<T>(state: AsyncState<T>, action: Action<T>): AsyncState<T> {
  switch (action.type) {
    case 'idle':
      // Returning the SAME object when already idle lets React bail out of the
      // re-render entirely, instead of looping.
      return state.status === 'idle' ? state : IDLE;
    case 'loading':
      return { status: 'loading', data: null, error: null };
    case 'success':
      return { status: 'success', data: action.data, error: null };
    case 'error':
      return { status: 'error', data: null, error: action.error };
  }
}

/**
 * Runs an async loader and reports its state, cancelling the in-flight request
 * whenever the dependencies change or the component unmounts.
 *
 * The AbortController is the important part. Without it, clicking 1914 then
 * 1969 quickly can leave the slower 1914 response arriving last and overwriting
 * the screen — a race condition that is invisible on a fast local network and
 * very visible on a phone. Aborting the previous request makes the bug
 * impossible rather than unlikely.
 *
 * `deps` is passed through to useEffect, so callers control when a reload happens.
 */
export function useAsyncResource<T>(
  load: (signal: AbortSignal) => Promise<T>,
  deps: readonly unknown[],
  options: { enabled?: boolean } = {},
): AsyncState<T> {
  const enabled = options.enabled ?? true;
  const [state, dispatch] = useReducer(reducer<T>, IDLE as AsyncState<T>);

  useEffect(() => {
    // A disabled resource must forget whatever it last loaded. Without this,
    // clearing the country selection would leave the PREVIOUS country's events
    // on screen — stale data that looks perfectly plausible to the user.
    if (!enabled) {
      dispatch({ type: 'idle' });
      return;
    }

    const controller = new AbortController();
    dispatch({ type: 'loading' });

    load(controller.signal)
      .then((data) => {
        if (!controller.signal.aborted) dispatch({ type: 'success', data });
      })
      .catch((error: unknown) => {
        if (controller.signal.aborted) return; // We cancelled it; not an error.
        dispatch({
          type: 'error',
          error:
            error instanceof ApiError
              ? error
              : new ApiError(0, 'UNKNOWN_ERROR', 'Something unexpected went wrong.'),
        });
      });

    return () => controller.abort();
    // `load` is intentionally excluded: callers pass an inline arrow function,
    // which is a new reference on every render and would loop forever. `deps`
    // is the explicit, honest list of what this request actually depends on.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [...deps, enabled]);

  return state;
}
