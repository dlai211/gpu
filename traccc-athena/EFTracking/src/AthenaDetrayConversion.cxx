
#include "AthenaDetrayConversion.h"
#include <sstream>
#include <stdexcept>

Identifier findAthenaID(
	  uint64_t detray_id
	, std::unordered_map<uint64_t, Identifier> const * detray_to_athena_map
	)
{
	if (detray_to_athena_map == nullptr) {
		throw std::runtime_error("detray to athena map null");
	}
	if (auto it = detray_to_athena_map->find(detray_id);
			it != detray_to_athena_map->end())
	{
		return it->second;
	} else {
		std::ostringstream oss;
		oss << "detray to athena lookup fail: " << detray_id;
		// throw std::runtime_error(oss.str());
		// std::cout << oss.str() << std::endl;
	}
	return Identifier();
}

uint64_t findDetrayID(
	  Identifier athena_id
	, std::unordered_map<Identifier, uint64_t> const * athena_to_detray_map
	)
{
	if (athena_to_detray_map == nullptr) {
		throw std::runtime_error("athena to detray map null");
	}
	if (auto it = athena_to_detray_map->find(athena_id);
			it != athena_to_detray_map->end())
	{
		return it->second;
	} else {
		std::ostringstream oss;
		oss << "detray to athena lookup fail: " << athena_id;
		// throw std::runtime_error(oss.str());
		// std::cout << oss.str() << std::endl;
	}
	return 0;
}