import type { Country, HistoricalEvent } from '../../domain/index.js';

/**
 * A tiny in-process "database" shared by the in-memory repositories.
 *
 * OOP lesson: COMPOSITION. Rather than each repository owning its own copy of
 * the data (and drifting apart), both repositories are *given* the same dataset
 * object. They compose it; they do not inherit from it. "Favour composition over
 * inheritance" in one small, concrete example.
 */
export class InMemoryDataset {
  constructor(
    readonly countries: readonly Country[] = [],
    readonly events: readonly HistoricalEvent[] = [],
  ) {}
}
