-- Migration 001: the two tables the MVP needs.
--
-- Design notes worth understanding:
--   * `code` is the natural key (ISO 3166-1 alpha-3) but we still keep a
--     surrogate integer `id`. Surrogate keys keep foreign keys narrow and let a
--     country be renamed or recoded without rewriting every child row.
--   * `year` is a plain INTEGER, not a DATE. Most historical events have no
--     reliable day, and DATE cannot represent "sometime in 1347". `month_day`
--     carries the extra precision on the occasions we have it.
--   * ON DELETE CASCADE: deleting a country deletes its events. Orphan events
--     would be unreachable rows that still inflate every count.

CREATE TABLE countries (
    id      SERIAL PRIMARY KEY,
    code    CHAR(3)      NOT NULL UNIQUE,
    name    TEXT         NOT NULL UNIQUE,
    region  TEXT         NOT NULL,
    CONSTRAINT countries_code_is_alpha CHECK (code ~ '^[A-Z]{3}$')
);

CREATE TABLE historical_events (
    id          SERIAL PRIMARY KEY,
    country_id  INTEGER NOT NULL REFERENCES countries (id) ON DELETE CASCADE,
    year        INTEGER NOT NULL,
    month_day   CHAR(5),
    title       TEXT    NOT NULL,
    summary     TEXT    NOT NULL,
    category    TEXT    NOT NULL,
    source_url  TEXT    NOT NULL,
    created_at  TIMESTAMPTZ NOT NULL DEFAULT now(),

    -- The database enforces the same invariants as the domain objects.
    -- Belt and braces: application code can be bypassed (psql, a bulk import,
    -- a future service in another language). The database cannot.
    CONSTRAINT events_year_in_range   CHECK (year BETWEEN 100 AND 2200),
    CONSTRAINT events_month_day_shape CHECK (month_day IS NULL OR month_day ~ '^\d{2}-\d{2}$'),
    CONSTRAINT events_title_length    CHECK (char_length(title) BETWEEN 1 AND 160),
    CONSTRAINT events_summary_length  CHECK (char_length(summary) BETWEEN 1 AND 1000),
    CONSTRAINT events_source_is_url   CHECK (source_url ~ '^https?://'),
    CONSTRAINT events_category_known  CHECK (
        category IN ('politics', 'conflict', 'science', 'culture', 'disaster', 'economy', 'society')
    ),
    -- The same event should not be recorded twice for a country in a year.
    CONSTRAINT events_unique_per_country_year UNIQUE (country_id, year, title)
);

-- The MVP's hot query is "events for (year, country)". A composite index in
-- this column order serves that AND the broader "events for (year)".
CREATE INDEX idx_events_year_country ON historical_events (year, country_id);
