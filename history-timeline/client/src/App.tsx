import { useCallback, useEffect, useState } from 'react';
import { Timeline, CountryPicker, EventList, Loading, ErrorMessage, EmptyMessage } from './components/index.js';
import { useCountriesForYear, useEvents, useTimeline } from './hooks/index.js';
import './styles/app.css';

/**
 * The app shell, and the single owner of "what is the user looking at".
 *
 * STATE DESIGN — the most important decision in a React app:
 * only two things are stored, `selectedYear` and `selectedCountryCode`. Every
 * other thing on screen is DERIVED from those two by a hook. There is no copy of
 * the country list, no `events` state to keep in sync, and therefore no way for
 * the screen to contradict itself.
 *
 * The rule of thumb: if you can compute it, do not store it.
 */
export function App() {
  const [selectedYear, setSelectedYear] = useState<number | null>(null);
  // Deriving the panels from `status` rather than from these two values is what
  // guarantees the screen can never show one year's events under another's label.
  const [selectedCountryCode, setSelectedCountryCode] = useState<string | null>(null);

  const timeline = useTimeline();
  const countries = useCountriesForYear(selectedYear);
  const events = useEvents(selectedYear, selectedCountryCode);

  // Start on the most recent year that actually has events, so a first-time
  // visitor sees content instead of an empty panel.
  useEffect(() => {
    if (selectedYear !== null || timeline.status !== 'success') return;
    const populated = [...timeline.data.ticks].reverse().find((t) => t.eventCount > 0);
    setSelectedYear(populated?.year ?? timeline.data.endYear);
  }, [timeline, selectedYear]);

  // Changing the year must clear the country: "France in 1969" is a valid pair,
  // but after switching to 1848 the same country may have nothing recorded.
  // Resetting dependent state when its parent changes is the cure for a whole
  // family of stale-UI bugs.
  const handleSelectYear = useCallback((year: number) => {
    setSelectedYear(year);
    setSelectedCountryCode(null);
  }, []);

  return (
    <div className="app">
      <header className="app__header">
        <h1 className="app__title">History Timeline</h1>
        <p className="app__tagline">
          Pick a year, choose a country, and read what is documented to have happened there. Every
          entry links to its source.
        </p>
      </header>

      <main className="app__main">
        <section className="panel panel--timeline" aria-labelledby="timeline-heading">
          <h2 id="timeline-heading" className="panel__heading">
            1. Choose a year
          </h2>
          {timeline.status === 'loading' && <Loading label="Loading the timeline…" />}
          {timeline.status === 'error' && <ErrorMessage error={timeline.error} />}
          {timeline.status === 'success' && (
            <Timeline
              ticks={timeline.data.ticks}
              selectedYear={selectedYear}
              onSelectYear={handleSelectYear}
            />
          )}
        </section>

        <div className="app__columns">
          <section className="panel panel--countries" aria-labelledby="countries-heading">
            <h2 id="countries-heading" className="panel__heading">
              2. Choose a country
              {selectedYear !== null && <span className="panel__context"> in {selectedYear}</span>}
            </h2>

            {countries.status === 'idle' && <EmptyMessage>Select a year to begin.</EmptyMessage>}
            {countries.status === 'loading' && <Loading label="Finding countries…" />}
            {countries.status === 'error' && <ErrorMessage error={countries.error} />}
            {countries.status === 'success' &&
              (countries.data.countries.length === 0 ? (
                <EmptyMessage>
                  Nothing is recorded for {selectedYear} yet. Try a year with a taller bar on the
                  timeline.
                </EmptyMessage>
              ) : (
                <CountryPicker
                  countries={countries.data.countries}
                  selectedCode={selectedCountryCode}
                  year={countries.data.year}
                  onSelectCountry={setSelectedCountryCode}
                />
              ))}
          </section>

          <section className="panel panel--events" aria-labelledby="events-heading">
            <h2 id="events-heading" className="panel__heading">
              3. What happened
            </h2>

            {events.status === 'idle' && (
              <EmptyMessage>Choose a country to see its events.</EmptyMessage>
            )}
            {events.status === 'loading' && <Loading label="Loading events…" />}
            {events.status === 'error' && <ErrorMessage error={events.error} />}
            {events.status === 'success' &&
              (events.data.events.length === 0 ? (
                <EmptyMessage>
                  Nothing is recorded for {events.data.country.name} in {events.data.year}.
                </EmptyMessage>
              ) : (
                <EventList
                  events={events.data.events}
                  countryName={events.data.country.name}
                  year={events.data.year}
                />
              ))}
          </section>
        </div>
      </main>

      <footer className="app__footer">
        <p>
          An educational project. Summaries are neutral descriptions of documented events; follow
          the source link on each card to read further.
        </p>
      </footer>
    </div>
  );
}
