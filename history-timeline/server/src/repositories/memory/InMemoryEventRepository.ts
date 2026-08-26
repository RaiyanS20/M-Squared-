import { EventRepository, type YearEventCount } from '../EventRepository.js';
import type { HistoricalEvent } from '../../domain/index.js';
import type { InMemoryDataset } from './InMemoryDataset.js';

export class InMemoryEventRepository extends EventRepository {
  constructor(private readonly data: InMemoryDataset) {
    super();
  }

  async findByYearAndCountry(year: number, countryId: number): Promise<HistoricalEvent[]> {
    return this.data.events
      .filter((e) => e.year.value === year && e.countryId === countryId)
      .sort((a, b) => a.compareByDate(b));
  }

  async findByYear(year: number): Promise<HistoricalEvent[]> {
    return this.data.events
      .filter((e) => e.year.value === year)
      .sort((a, b) => a.compareByDate(b));
  }

  async countByYearRange(startYear: number, endYear: number): Promise<YearEventCount[]> {
    const counts = new Map<number, number>();
    for (const event of this.data.events) {
      const y = event.year.value;
      if (y >= startYear && y <= endYear) {
        counts.set(y, (counts.get(y) ?? 0) + 1);
      }
    }
    return [...counts.entries()]
      .map(([year, eventCount]) => ({ year, eventCount }))
      .sort((a, b) => a.year - b.year);
  }
}
