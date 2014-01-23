#ifndef FIBERIO_H
#define FIBERIO_H

#include "dtitypes.h"

GroupType::Pointer readFiberFile(const std::string & filename);

void writeFiberFile(const std::string & filename, GroupType::Pointer fibergroup, bool saveProperties = true);

#endif
