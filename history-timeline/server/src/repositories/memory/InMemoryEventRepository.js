import { EventRepository } from '../EventRepository.js';

export class InMemoryEventRepository extends EventRepository {
  #data;

  constructor(dataset) {
    super();
    this.#data = dataset;
  }

  async findByYearAndCountry(year, countryId) {
    return this.#data.events
      .filter((e) => e.year.value === year && e.countryId === countryId)
      .sort((a, b) => a.compareByDate(b));
  }

  async findByYear(year) {
    return this.#data.events.filter((e) => e.year.value === year).sort((a, b) => a.compareByDate(b));
  }

  async countByYearRange(startYear, endYear) {
    const counts = new Map();
    for (const event of this.#data.events) {
      const y = event.year.value;
      if (y >= startYear && y <= endYear) counts.set(y, (counts.get(y) ?? 0) + 1);
    }
    return [...counts.entries()]
      .map(([year, eventCount]) => ({ year, eventCount }))
      .sort((a, b) => a.year - b.year);
  }
}
