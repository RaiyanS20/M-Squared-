import { describe, expect, it, vi } from 'vitest';
import { ApiError, TimelineApi } from '../src/api/index.js';
import { COUNTRIES_1969, TIMELINE, stubFetch } from './stubApi.js';

describe('TimelineApi', () => {
  it('requests the timeline and returns the parsed body', async () => {
    const fetchSpy = stubFetch({ '/api/timeline': TIMELINE });
    const api = new TimelineApi('');

    await expect(api.getTimeline()).resolves.toEqual(TIMELINE);
    expect(fetchSpy).toHaveBeenCalledWith('/api/timeline', expect.anything());
  });

  it('appends only the query parameters that were provided', async () => {
    const fetchSpy = stubFetch({ '/api/timeline': TIMELINE });
    await new TimelineApi('').getTimeline(1900, undefined);
    expect(fetchSpy.mock.calls[0]![0]).toBe('/api/timeline?startYear=1900');
  });

  it('unwraps the countries array', async () => {
    stubFetch({ '/api/countries': { countries: COUNTRIES_1969.countries } });
    await expect(new TimelineApi('').getAllCountries()).resolves.toHaveLength(2);
  });

  it('honours a configured base URL', async () => {
    const fetchSpy = stubFetch({ '/api/timeline': TIMELINE });
    await new TimelineApi('https://api.example.com').getTimeline();
    expect(fetchSpy.mock.calls[0]![0]).toBe('https://api.example.com/api/timeline');
  });

  describe('error handling', () => {
    it('turns a 404 body into an ApiError carrying the server code', async () => {
      stubFetch(
        { '/api/years': { error: { code: 'NOT_FOUND', message: "Country 'ZZZ' was not found." } } },
        { status: { '/api/years': 404 } },
      );

      const error = await new TimelineApi('').getEvents(1969, 'ZZZ').catch((e: unknown) => e);
      expect(error).toBeInstanceOf(ApiError);
      expect((error as ApiError).code).toBe('NOT_FOUND');
      expect((error as ApiError).status).toBe(404);
      expect((error as ApiError).isClientError).toBe(true);
    });

    it('survives an error response that is not JSON', async () => {
      vi.spyOn(globalThis, 'fetch').mockResolvedValue(
        new Response('<html>502 Bad Gateway</html>', { status: 502 }),
      );

      const error = (await new TimelineApi('').getTimeline().catch((e) => e)) as ApiError;
      expect(error).toBeInstanceOf(ApiError);
      expect(error.code).toBe('UNKNOWN_ERROR');
      expect(error.isClientError).toBe(false); // 5xx — worth retrying.
    });

    it('reports a network failure as NETWORK_ERROR', async () => {
      vi.spyOn(globalThis, 'fetch').mockRejectedValue(new TypeError('Failed to fetch'));

      const error = (await new TimelineApi('').getTimeline().catch((e) => e)) as ApiError;
      expect(error.code).toBe('NETWORK_ERROR');
      expect(error.message).toMatch(/could not reach/i);
    });

    it('re-throws an abort untouched so callers can ignore it', async () => {
      vi.spyOn(globalThis, 'fetch').mockRejectedValue(
        new DOMException('The operation was aborted.', 'AbortError'),
      );

      const error = await new TimelineApi('').getTimeline().catch((e: unknown) => e);
      expect(error).toBeInstanceOf(DOMException);
      expect(error).not.toBeInstanceOf(ApiError);
    });
  });
});
