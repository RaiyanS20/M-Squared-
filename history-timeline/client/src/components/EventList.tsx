import type { HistoricalEvent } from '../api/index.js';

/** Turns "04-27" into "27 April" for display. Pure, so it is trivially testable. */
export function formatMonthDay(monthDay: string | null): string | null {
  if (!monthDay) return null;
  const [month, day] = monthDay.split('-').map(Number);
  if (!month || !day || month < 1 || month > 12) return null;
  const monthName = new Date(Date.UTC(2000, month - 1, 1)).toLocaleString('en-GB', {
    month: 'long',
    timeZone: 'UTC',
  });
  return `${day} ${monthName}`;
}

interface EventListProps {
  events: HistoricalEvent[];
  countryName: string;
  year: number;
}

/**
 * MVP feature 4: the payoff.
 *
 * Note the source link on every card. That is not decoration — it is the whole
 * editorial promise of the project made visible, and the domain layer refuses to
 * create an event that could not render one.
 */
export function EventList({ events, countryName, year }: EventListProps) {
  return (
    <ol className="events" aria-label={`Events in ${countryName} in ${year}`}>
      {events.map((event) => {
        const date = formatMonthDay(event.monthDay);
        return (
          <li key={event.id} className="events__item">
            <article className="event">
              <header className="event__header">
                <span className={`event__category event__category--${event.category}`}>
                  {event.category}
                </span>
                <time className="event__date" dateTime={dateTimeAttr(event.year, event.monthDay)}>
                  {date ? `${date} ${event.year}` : event.year}
                </time>
              </header>
              <h3 className="event__title">{event.title}</h3>
              <p className="event__summary">{event.summary}</p>
              <a
                className="event__source"
                href={event.sourceUrl}
                target="_blank"
                // `noopener` stops the opened page reaching back through
                // window.opener; `noreferrer` withholds the referring URL.
                rel="noopener noreferrer"
              >
                Read the source
                <span className="visually-hidden"> for {event.title} (opens in a new tab)</span>
              </a>
            </article>
          </li>
        );
      })}
    </ol>
  );
}

/** A machine-readable value for <time>: "1994-04-27", or just "1994". */
function dateTimeAttr(year: number, monthDay: string | null): string {
  return monthDay ? `${year}-${monthDay}` : String(year);
}
