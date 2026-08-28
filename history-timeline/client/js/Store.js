/**
 * A tiny observable state container — about 40 lines.
 *
 * This is the piece a framework would give you (React's `useState`, Vue's
 * reactivity, Redux). Written out, the whole idea is:
 *
 *   1. hold some state,
 *   2. let interested parties SUBSCRIBE to changes,
 *   3. when state changes, NOTIFY them.
 *
 * That is the observer pattern, and once you have seen it in 40 lines, every
 * state library you meet later is a variation on it.
 *
 * WHY BOTHER, instead of just setting variables and updating the DOM inline?
 * Because with scattered variables there is no single answer to "what is on
 * screen right now?", and two parts of the UI eventually disagree — the heading
 * says 1994 while the list below still shows 1969. One state object, one place
 * that changes it, and every view derived from it, makes that impossible.
 */
export class Store {
  #state;
  #listeners = new Set();
  #notifying = false;

  constructor(initialState = {}) {
    this.#state = Object.freeze({ ...initialState });
  }

  /** The current state. Frozen, so a view cannot mutate it by accident. */
  getState() {
    return this.#state;
  }

  /**
   * Merges a patch into state and notifies subscribers.
   *
   * State is REPLACED, never mutated: `{ ...old, ...patch }` builds a new
   * object. That means a subscriber can compare `previous.events !== next.events`
   * with `!==` to know whether that slice actually changed — which is exactly
   * how the views below avoid rebuilding the whole page on every keystroke.
   */
  setState(patch) {
    const next = Object.freeze({ ...this.#state, ...patch });

    // Nothing to do if every value is identical. Cheap, and it stops a
    // subscriber that calls setState from looping forever.
    const unchanged = Object.keys(patch).every((key) => this.#state[key] === next[key]);
    if (unchanged) return;

    const previous = this.#state;
    this.#state = next;

    if (this.#notifying) {
      // Setting state from inside a listener is a design smell and, unguarded,
      // an infinite loop. Fail loudly rather than hanging the browser tab.
      throw new Error('setState was called while notifying subscribers.');
    }

    this.#notifying = true;
    try {
      // Copy the set before iterating: a listener that unsubscribes during
      // notification would otherwise mutate the collection being iterated.
      for (const listener of [...this.#listeners]) listener(next, previous);
    } finally {
      this.#notifying = false;
    }
  }

  /** Subscribes to changes. Returns a function that unsubscribes. */
  subscribe(listener) {
    this.#listeners.add(listener);
    return () => this.#listeners.delete(listener);
  }
}

/**
 * The four states any remote request can be in.
 *
 * Compare the usual alternative — three separate variables, `data`, `loading`
 * and `error` — which allows `loading === true` AND `error` set AND stale `data`
 * all at once. Eight combinations, half of them nonsense.
 *
 * One `status` field permits only the four states that are real, and reading
 * `state.events.status` tells you everything.
 */
export const Async = Object.freeze({
  /** Not requested — usually because a prerequisite has not been chosen yet. */
  idle: () => ({ status: 'idle', data: null, error: null }),
  loading: () => ({ status: 'loading', data: null, error: null }),
  success: (data) => ({ status: 'success', data, error: null }),
  error: (error) => ({ status: 'error', data: null, error }),
});
