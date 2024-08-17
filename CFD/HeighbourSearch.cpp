#include "pch.h"
#include "HeighbourSearch.h"

uint32_t getHash(const SPH::int2& cell)
{
    return (
        (unsigned int)(cell.x * 73856093)
        ^ (unsigned int)(cell.y * 19349663)
        ) % TABLE_SIZE;
}

SPH::int2 getCell(std::shared_ptr<SPH::SParticle> p, float h)
{
    return SPH::int2(p->r().x / h, p->r().y / h);
}

uint32_t* createNeighborTable(std::vector<std::shared_ptr<SPH::SParticle>>& particles)
{
    uint32_t* particleTable
        = (uint32_t*)malloc(sizeof(uint32_t) * TABLE_SIZE);
    for (size_t i = 0; i < TABLE_SIZE; ++i) {
        particleTable[i] = NO_PARTICLE;
    }

    uint32_t prevHash = NO_PARTICLE;
    for (size_t i = 0; i < particles.size(); ++i) {
        uint32_t currentHash = particles[i]->hash;
        if (currentHash != prevHash) {
            particleTable[currentHash] = i;
            prevHash = currentHash;
        }
    }
    return particleTable;
}
