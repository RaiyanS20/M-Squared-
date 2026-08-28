/**
 * A tiny in-process "database" shared by the in-memory repositories.
 *
 * OOP lesson: COMPOSITION. Rather than each repository owning its own copy of
 * the data (and drifting apart), both are *given* the same dataset object. They
 * compose it; they do not inherit from it.
 *
 * "Favour composition over inheritance" in one concrete example: the
 * relationship here is "a repository USES a dataset", not "a repository IS a
 * dataset" — so it is a constructor argument, not a base class.
 */
export class InMemoryDataset {
  constructor(countries = [], events = []) {
    this.countries = countries;
    this.events = events;
    Object.freeze(this);
  }
}
