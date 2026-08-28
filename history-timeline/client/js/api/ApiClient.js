/**
 * An error carrying the server's machine-readable code, so the UI can react to
 * WHAT went wrong rather than pattern-matching on English message text.
 *
 * Extending the built-in `Error` means `instanceof`, stack traces and
 * `console.error` all behave normally.
 */
export class ApiError extends Error {
  constructor(status, code, message) {
    super(message);
    this.name = 'ApiError';
    this.status = status;
    this.code = code;
  }

  /** 4xx means the request was wrong; retrying it unchanged will not help. */
  get isClientError() {
    return this.status >= 400 && this.status < 500;
  }
}

/**
 * Base class for anything that talks HTTP to our backend.
 *
 * OOP lesson: this is where classes clearly earn their place on the front end.
 * The network layer has state (a base URL), behaviour (request, parse, raise)
 * and a natural subclass per API area. Inheriting from `ApiClient` means
 * `TimelineApi` never repeats error handling or JSON parsing.
 *
 * No other file in the browser app calls `fetch`. Adding authentication,
 * retries or a different base URL is therefore a change to ONE file — the same
 * principle as services never writing SQL on the server.
 */
export class ApiClient {
  #baseUrl;

  constructor(baseUrl = '') {
    if (new.target === ApiClient) {
      throw new TypeError('ApiClient is abstract; extend it.');
    }
    this.#baseUrl = baseUrl;
  }

  /**
   * One place where every HTTP concern is handled: URL building, headers,
   * non-2xx responses, malformed JSON, and network failure.
   *
   * A detail people miss: `fetch` only rejects on NETWORK failure. A 404 or a
   * 500 is a perfectly successful fetch with `response.ok === false`. Forgetting
   * to check `.ok` is the most common bug in front-end HTTP code.
   */
  async get(path, signal) {
    let response;
    try {
      response = await fetch(`${this.#baseUrl}${path}`, {
        headers: { Accept: 'application/json' },
        signal,
      });
    } catch (cause) {
      // An aborted request is a normal part of the lifecycle (the user moved
      // on), not a failure to report. Re-throw it untouched so callers can
      // recognise and ignore it.
      if (cause?.name === 'AbortError') throw cause;
      throw new ApiError(0, 'NETWORK_ERROR', 'Could not reach the history service.');
    }

    if (!response.ok) {
      // Our API always returns { error: { code, message } } — but a proxy or a
      // crashed process might return HTML, so parsing must not itself throw.
      const body = await response.json().catch(() => null);
      throw new ApiError(
        response.status,
        body?.error?.code ?? 'UNKNOWN_ERROR',
        body?.error?.message ?? `Request failed with status ${response.status}.`,
      );
    }

    return response.json();
  }

  /** Builds a query string, omitting undefined values. */
  query(params) {
    const search = new URLSearchParams();
    for (const [key, value] of Object.entries(params)) {
      if (value !== undefined && value !== null) search.set(key, String(value));
    }
    const s = search.toString();
    return s ? `?${s}` : '';
  }
}
