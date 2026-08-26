import { useMemo } from 'react';
import { timelineApi } from '../api/index.js';
import type { CountriesForYear, EventsForYearAndCountry, TimelineScale } from '../api/index.js';
import { useAsyncResource, type AsyncState } from './useAsyncResource.js';

/**
 * One hook per MVP feature. Each is a thin, named wrapper around
 * `useAsyncResource` — components ask for "the countries in 1969", not for a
 * URL, and never see `fetch`, loading flags, or abort controllers.
 */

/** Feature 1. */
export function useTimeline(startYear?: number, endYear?: number): AsyncState<TimelineScale> {
  return useAsyncResource(
    (signal) => timelineApi.getTimeline(startYear, endYear, signal),
    [startYear, endYear],
  );
}

/** Features 2 + 3: only runs once a year has been chosen. */
export function useCountriesForYear(year: number | null): AsyncState<CountriesForYear> {
  return useAsyncResource(
    (signal) => timelineApi.getCountriesForYear(year!, signal),
    [year],
    { enabled: year !== null },
  );
}

/** Feature 4: only runs once both a year and a country have been chosen. */
export function useEvents(
  year: number | null,
  countryCode: string | null,
): AsyncState<EventsForYearAndCountry> {
  return useAsyncResource(
    (signal) => timelineApi.getEvents(year!, countryCode!, signal),
    [year, countryCode],
    { enabled: year !== null && countryCode !== null },
  );
}

/** Groups countries by region so the list is browsable rather than a wall of names. */
export function useCountriesByRegion<T extends { region: string; name: string }>(
  countries: readonly T[] | undefined,
): Array<[string, T[]]> {
  return useMemo(() => {
    if (!countries) return [];
    const byRegion = new Map<string, T[]>();
    for (const country of countries) {
      const bucket = byRegion.get(country.region) ?? [];
      bucket.push(country);
      byRegion.set(country.region, bucket);
    }
    // Sort regions, and countries within each region, so the list is stable and
    // scannable no matter what order the API happened to return.
    return [...byRegion.entries()]
      .map(([region, list]): [string, T[]] => [
        region,
        [...list].sort((a, b) => a.name.localeCompare(b.name)),
      ])
      .sort(([a], [b]) => a.localeCompare(b));
  }, [countries]);
}
