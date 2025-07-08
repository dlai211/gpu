
#include "Identifier/Identifier.h"

Identifier findAthenaID(
	  uint64_t detray_id
	, std::unordered_map<uint64_t, Identifier> const * detray_to_athena_map
	);

uint64_t findDetrayID(
	  Identifier athena_id
	, std::unordered_map<Identifier, uint64_t> const * athena_to_detray_map
	);