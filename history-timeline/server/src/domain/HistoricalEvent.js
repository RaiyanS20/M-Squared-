import { ValidationError } from './errors.js';
import { Year } from './Year.js';

/**
 * How we categorise an event.
 *
 * JavaScript has no `enum`, so this is the idiomatic stand-in: a frozen object
 * of constants. `Object.freeze` means a typo like `EventCategory.Politcs`
 * evaluates to `undefined` (and fails validation below) rather than silently
 * adding a new category.
 */
export const EventCategory = Object.freeze({
  Politics: 'politics',
  Conflict: 'conflict',
  Science: 'science',
  Culture: 'culture',
  Disaster: 'disaster',
  Economy: 'economy',
  Society: 'society',
});

/** Every valid category value, for validation and for the UI to iterate. */
export const EVENT_CATEGORIES = Object.freeze(Object.values(EventCategory));

/**
 * A single factual event that happened in one country in one year.
 *
 * OOP lesson: an entity is not a bag of data — it is data PLUS the rules that
 * keep that data honest. Every invariant this project cares about lives here, in
 * one constructor, instead of being re-checked in the controller, the seed
 * script and the browser.
 */
export class HistoricalEvent {
  static MAX_TITLE_LENGTH = 160;
  static MAX_SUMMARY_LENGTH = 1000;

  constructor({ id, countryId, year, title, summary, category, sourceUrl, monthDay = null }) {
    if (!(year instanceof Year)) {
      // Accepting a raw number here would let an unvalidated year in through the
      // back door, defeating the point of having a Year class at all.
      throw new ValidationError('Event year must be a Year instance.', 'year');
    }

    const cleanTitle = typeof title === 'string' ? title.trim() : '';
    if (!cleanTitle) {
      throw new ValidationError('Event title is required.', 'title');
    }
    if (cleanTitle.length > HistoricalEvent.MAX_TITLE_LENGTH) {
      throw new ValidationError(
        `Event title must be at most ${HistoricalEvent.MAX_TITLE_LENGTH} characters.`,
        'title',
      );
    }

    const cleanSummary = typeof summary === 'string' ? summary.trim() : '';
    if (!cleanSummary) {
      throw new ValidationError('Event summary is required.', 'summary');
    }
    if (cleanSummary.length > HistoricalEvent.MAX_SUMMARY_LENGTH) {
      throw new ValidationError(
        `Event summary must be at most ${HistoricalEvent.MAX_SUMMARY_LENGTH} characters.`,
        'summary',
      );
    }

    if (!EVENT_CATEGORIES.includes(category)) {
      throw new ValidationError(`Unknown event category ${JSON.stringify(category)}.`, 'category');
    }

    // THE EDITORIAL INVARIANT OF THE WHOLE PROJECT.
    // This is a factual, educational history site, so an event without a
    // citation is a rumour — and the domain refuses to construct one. Not
    // "should not": cannot. Not from the API, not from the seed, not in a test.
    if (typeof sourceUrl !== 'string' || !/^https?:\/\/\S+$/.test(sourceUrl)) {
      throw new ValidationError(
        `Event must cite a source URL, received ${JSON.stringify(sourceUrl)}.`,
        'sourceUrl',
      );
    }

    if (monthDay != null && !/^\d{2}-\d{2}$/.test(monthDay)) {
      throw new ValidationError(
        `monthDay must look like "MM-DD", received ${JSON.stringify(monthDay)}.`,
        'monthDay',
      );
    }

    this.id = id;
    this.countryId = countryId;
    this.year = year;
    this.title = cleanTitle;
    this.summary = cleanSummary;
    this.category = category;
    this.sourceUrl = sourceUrl;
    this.monthDay = monthDay ?? null;
    Object.freeze(this);
  }

  equals(other) {
    return other instanceof HistoricalEvent && other.id === this.id;
  }

  /**
   * Sort order: by year, then by date, and events with a known date come before
   * events that only have a year. Most of history is undated, so this is a real
   * requirement rather than an edge case.
   */
  compareByDate(other) {
    if (this.year.value !== other.year.value) return this.year.value - other.year.value;
    if (this.monthDay && other.monthDay) return this.monthDay.localeCompare(other.monthDay);
    if (this.monthDay) return -1;
    if (other.monthDay) return 1;
    return this.title.localeCompare(other.title);
  }

  toJSON() {
    return {
      id: this.id,
      countryId: this.countryId,
      year: this.year.value,
      title: this.title,
      summary: this.summary,
      category: this.category,
      sourceUrl: this.sourceUrl,
      monthDay: this.monthDay,
    };
  }
}
