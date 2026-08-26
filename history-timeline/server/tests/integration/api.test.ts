import { describe, expect, it } from 'vitest';
import request from 'supertest';
import type { Express } from 'express';
import { createApp } from '../../src/api/index.js';
import {
  InMemoryCountryRepository,
  InMemoryEventRepository,
} from '../../src/repositories/memory/index.js';
import { TimelineService } from '../../src/services/index.js';
import { aDataset } from '../fixtures.js';
import { Year } from '../../src/domain/index.js';

/**
 * INTEGRATION TESTS: the full HTTP stack — routing, validation, controllers,
 * services, repositories, error handling — with in-memory data underneath.
 *
 * This is the layer of the testing pyramid that catches the mistakes unit tests
 * cannot see: a route registered at the wrong path, a 500 where a 404 belongs,
 * a field renamed in the JSON that the React app still expects.
 *
 * Supertest binds an ephemeral port for the duration of each request, so these
 * still need no running server and no database.
 */
function buildTestApp(): Express {
  const data = aDataset();
  const service = new TimelineService(
    new InMemoryCountryRepository(data),
    new InMemoryEventRepository(data),
  );
  return createApp({ service, db: null, corsOrigins: ['http://localhost:5173'], isProduction: false });
}

describe('HTTP API', () => {
  const app = buildTestApp();

  describe('GET /health and /ready', () => {
    it('reports ok', async () => {
      const res = await request(app).get('/health').expect(200);
      expect(res.body.status).toBe('ok');
    });

    it('reports ready when no database is configured', async () => {
      const res = await request(app).get('/ready').expect(200);
      expect(res.body).toMatchObject({ status: 'ready', database: 'not-configured' });
    });
  });

  describe('GET /api/timeline', () => {
    it('returns the default MVP scale', async () => {
      const res = await request(app).get('/api/timeline').expect(200);
      expect(res.body.startYear).toBe(1900);
      expect(res.body.endYear).toBe(Year.latestAllowed());
      expect(res.body.ticks[0]).toEqual({ year: 1900, eventCount: 0 });
    });

    it('honours an explicit range', async () => {
      const res = await request(app).get('/api/timeline?startYear=1969&endYear=1970').expect(200);
      expect(res.body.ticks).toEqual([
        { year: 1969, eventCount: 3 },
        { year: 1970, eventCount: 0 },
      ]);
    });

    it('rejects a non-numeric year with 400 and a machine-readable code', async () => {
      const res = await request(app).get('/api/timeline?startYear=abc').expect(400);
      expect(res.body.error.code).toBe('VALIDATION_ERROR');
      expect(res.body.error.message).toMatch(/not a valid year/i);
    });

    it('rejects a backwards range with 400', async () => {
      await request(app).get('/api/timeline?startYear=2000&endYear=1900').expect(400);
    });
  });

  describe('GET /api/countries', () => {
    it('lists every country', async () => {
      const res = await request(app).get('/api/countries').expect(200);
      expect(res.body.countries.map((c: { code: string }) => c.code)).toEqual(['BRA', 'FRA', 'JPN']);
    });

    it('exposes only the intended fields', async () => {
      const res = await request(app).get('/api/countries').expect(200);
      expect(Object.keys(res.body.countries[0]).sort()).toEqual(['code', 'id', 'name', 'region']);
    });
  });

  describe('GET /api/years/:year/countries', () => {
    it('lists only countries with events in that year', async () => {
      const res = await request(app).get('/api/years/1969/countries').expect(200);
      expect(res.body.year).toBe(1969);
      expect(res.body.countries.map((c: { code: string }) => c.code)).toEqual(['FRA', 'JPN']);
    });

    it('returns an empty list for a quiet year', async () => {
      const res = await request(app).get('/api/years/1902/countries').expect(200);
      expect(res.body.countries).toEqual([]);
    });

    it('rejects a year before 100 AD', async () => {
      await request(app).get('/api/years/50/countries').expect(400);
    });

    it('rejects a future year', async () => {
      await request(app).get(`/api/years/${Year.latestAllowed() + 1}/countries`).expect(400);
    });
  });

  describe('GET /api/years/:year/countries/:code/events', () => {
    it('returns the events for that year and country, in date order', async () => {
      const res = await request(app).get('/api/years/1969/countries/FRA/events').expect(200);
      expect(res.body.year).toBe(1969);
      expect(res.body.country.name).toBe('France');
      expect(res.body.events.map((e: { title: string }) => e.title)).toEqual([
        'Concorde first flight',
        'De Gaulle resigns',
      ]);
    });

    it('serialises an event with the fields the client needs', async () => {
      const res = await request(app).get('/api/years/1969/countries/FRA/events').expect(200);
      expect(Object.keys(res.body.events[0]).sort()).toEqual([
        'category', 'countryId', 'id', 'monthDay', 'sourceUrl', 'summary', 'title', 'year',
      ]);
    });

    it('accepts a lowercase code', async () => {
      await request(app).get('/api/years/1969/countries/fra/events').expect(200);
    });

    it('returns 404 for an unknown country', async () => {
      const res = await request(app).get('/api/years/1969/countries/ZZZ/events').expect(404);
      expect(res.body.error.code).toBe('NOT_FOUND');
    });

    it('returns 400 for a malformed country code', async () => {
      const res = await request(app).get('/api/years/1969/countries/FRANCE/events').expect(400);
      expect(res.body.error.code).toBe('VALIDATION_ERROR');
    });

    it('returns 200 with an empty list when nothing was recorded', async () => {
      const res = await request(app).get('/api/years/1955/countries/FRA/events').expect(200);
      expect(res.body.events).toEqual([]);
    });
  });

  describe('unknown routes', () => {
    it('return JSON, not an HTML error page', async () => {
      const res = await request(app).get('/api/nope').expect(404);
      expect(res.body.error.code).toBe('ROUTE_NOT_FOUND');
    });
  });

  describe('security headers and CORS', () => {
    it('does not advertise the framework', async () => {
      const res = await request(app).get('/health');
      expect(res.headers['x-powered-by']).toBeUndefined();
    });

    it('allows the configured origin', async () => {
      const res = await request(app).get('/api/countries').set('Origin', 'http://localhost:5173');
      expect(res.headers['access-control-allow-origin']).toBe('http://localhost:5173');
    });

    it('does not allow an unconfigured origin', async () => {
      const res = await request(app).get('/api/countries').set('Origin', 'https://evil.example');
      expect(res.headers['access-control-allow-origin']).toBeUndefined();
    });
  });
});
