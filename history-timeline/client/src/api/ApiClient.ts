/**
 * An error carrying the server's machine-readable code, so the UI can react to
 * *what* went wrong rather than pattern-matching on English message text.
 */
export class ApiError extends Error {
  constructor(
    readonly status: number,
    readonly code: string,
    message: string,
  ) {
    super(message);
    this.name = 'ApiError';
  }

  /** A 4xx means the request was wrong; retrying it unchanged will not help. */
  get isClientError(): boolean {
    return this.status >= 400 && this.status < 500;
  }
}

/**
 * Base class for anything that talks HTTP to our backend.
 *
 * OOP lesson: this is where classes still earn their place in modern React.
 * Components are functions — that is idiomatic React and we follow it — but the
 * NETWORK LAYER is a natural object: it has state (base URL, an abort signal),
 * behaviour (request, parse, raise), and a subclass per API area. Inheriting
 * from `ApiClient` means `TimelineApi` never repeats error handling or JSON
 * parsing, and a future `AdminApi` gets the same behaviour for free.
 */
export abstract class ApiClient {
  protected constructor(private readonly baseUrl: string) {}

  /**
   * One place where every HTTP concern is handled: URL building, headers,
   * non-2xx responses, malformed JSON, and network failure. Everything a
   * subclass writes is a one-line, fully typed method.
   */
  protected async get<T>(path: string, signal?: AbortSignal): Promise<T> {
    let response: Response;
    try {
      response = await fetch(`${this.baseUrl}${path}`, {
        headers: { Accept: 'application/json' },
        signal,
      });
    } catch (cause) {
      // An aborted request is a normal part of React's lifecycle (the user moved
      // on), not a failure to report. Re-throw it untouched so callers can ignore it.
      if (cause instanceof DOMException && cause.name === 'AbortError') throw cause;
      throw new ApiError(0, 'NETWORK_ERROR', 'Could not reach the history service.');
    }

    if (!response.ok) {
      // The API always returns { error: { code, message } } — but a proxy or a
      // crashed process might return HTML, so parsing must not itself throw.
      const body = (await response.json().catch(() => null)) as
        | { error?: { code?: string; message?: string } }
        | null;
      throw new ApiError(
        response.status,
        body?.error?.code ?? 'UNKNOWN_ERROR',
        body?.error?.message ?? `Request failed with status ${response.status}.`,
      );
    }

    return (await response.json()) as T;
  }

  /** Builds a query string, omitting undefined values. */
  protected query(params: Record<string, string | number | undefined>): string {
    const search = new URLSearchParams();
    for (const [key, value] of Object.entries(params)) {
      if (value !== undefined) search.set(key, String(value));
    }
    const s = search.toString();
    return s ? `?${s}` : '';
  }
}
