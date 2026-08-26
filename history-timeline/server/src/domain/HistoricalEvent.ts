import { ValidationError } from './errors.js';
import { Year } from './Year.js';

/**
 * How we categorise an event. A string enum (rather than a bare string) means
 * the compiler catches `'Politcal'` typos, and the UI can exhaustively switch
 * over every case.
 */
export enum EventCategory {
  Politics = 'politics',
  Conflict = 'conflict',
  Science = 'science',
  Culture = 'culture',
  Disaster = 'disaster',
  Economy = 'economy',
  Society = 'society',
}

export interface HistoricalEventProps {
  id: number;
  countryId: number;
  year: Year;
  title: string;
  summary: string;
  category: EventCategory;
  /** Where a curious learner can verify the claim. Required by our editorial rule. */
  sourceUrl: string;
  /** Optional day-level precision, ISO "MM-DD". Many old events only have a year. */
  monthDay?: string | null;
}

/**
 * A single factual event that happened in one country in one year.
 *
 * OOP lesson: an entity is not a bag of data — it is data PLUS the rules that
 * keep that data honest. Every invariant this project cares about ("an event
 * must cite a source", "the summary must be short enough to skim") lives here,
 * in one place, instead of being re-checked in the controller, the seed script
 * and the React form.
 */
export class HistoricalEvent {
  static readonly MAX_TITLE_LENGTH = 160;
  static readonly MAX_SUMMARY_LENGTH = 1000;

  readonly id: number;
  readonly countryId: number;
  readonly year: Year;
  readonly title: string;
  readonly summary: string;
  readonly category: EventCategory;
  readonly sourceUrl: string;
  readonly monthDay: string | null;

  constructor(props: HistoricalEventProps) {
    const title = props.title?.trim() ?? '';
    if (!title) {
      throw new ValidationError('Event title is required.', 'title');
    }
    if (title.length > HistoricalEvent.MAX_TITLE_LENGTH) {
      throw new ValidationError(
        `Event title must be at most ${HistoricalEvent.MAX_TITLE_LENGTH} characters.`,
        'title',
      );
    }

    const summary = props.summary?.trim() ?? '';
    if (!summary) {
      throw new ValidationError('Event summary is required.', 'summary');
    }
    if (summary.length > HistoricalEvent.MAX_SUMMARY_LENGTH) {
      throw new ValidationError(
        `Event summary must be at most ${HistoricalEvent.MAX_SUMMARY_LENGTH} characters.`,
        'summary',
      );
    }

    if (!Object.values(EventCategory).includes(props.category)) {
      throw new ValidationError(`Unknown event category "${props.category}".`, 'category');
    }

    // Editorial invariant: this is an educational, factual project. An event
    // without a citation is a rumour, so the domain refuses to construct one.
    if (!/^https?:\/\/\S+$/.test(props.sourceUrl ?? '')) {
      throw new ValidationError(
        `Event must cite a source URL, received "${props.sourceUrl}".`,
        'sourceUrl',
      );
    }

    if (props.monthDay != null && !/^\d{2}-\d{2}$/.test(props.monthDay)) {
      throw new ValidationError(
        `monthDay must look like "MM-DD", received "${props.monthDay}".`,
        'monthDay',
      );
    }

    this.id = props.id;
    this.countryId = props.countryId;
    this.year = props.year;
    this.title = title;
    this.summary = summary;
    this.category = props.category;
    this.sourceUrl = props.sourceUrl;
    this.monthDay = props.monthDay ?? null;
    Object.freeze(this);
  }

  equals(other: HistoricalEvent): boolean {
    return other instanceof HistoricalEvent && other.id === this.id;
  }

  /** Events with a known date sort before events that only have a year. */
  compareByDate(other: HistoricalEvent): number {
    if (this.year.value !== other.year.value) {
      return this.year.value - other.year.value;
    }
    if (this.monthDay && other.monthDay) {
      return this.monthDay.localeCompare(other.monthDay);
    }
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
