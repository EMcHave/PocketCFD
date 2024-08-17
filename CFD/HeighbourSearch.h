#pragma once

#include "SPH/SParticle.h"

const uint32_t TABLE_SIZE = 262144;
const uint32_t NO_PARTICLE = 0xFFFFFFFF;

/// Returns a hash of the cell position
uint32_t getHash(const SPH::int2& cell);

/// Get the cell that the particle is in.
SPH::int2 getCell(std::shared_ptr<SPH::SParticle> p, float h);

/// Creates the particle neighbor hash table.
/// It is the caller's responsibility to free the table.
uint32_t* createNeighborTable(
    std::vector<std::shared_ptr<SPH::SParticle> >& particles);
