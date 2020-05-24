#pragma once

#include <functional>

struct Data;
struct HeuristicsResults;

namespace nmHeuristics
{
	void OTR(const Data &data, const bool bRandom, HeuristicsResults &res);
	void OTRall(const Data &data, HeuristicsResults &res);
	void ORT(const Data &data, const bool bRandom, HeuristicsResults &res);
	void TRO(const Data &data, const bool bRandom, HeuristicsResults &res);
	void TOR(const Data &data, const bool bRandom, HeuristicsResults &res);
	void ROT(const Data &data, const bool bRandom, HeuristicsResults &res);
	void RTO(const Data &data, const bool bRandom, HeuristicsResults &res);
}

struct SHeuristic
{
	SHeuristic(std::function<void(const Data&, const bool, HeuristicsResults&)> fHeuristic, std::string sName)
		: fHeuristic(fHeuristic), sName(sName) { }

	const std::function<void(const Data&, const bool, HeuristicsResults&)> fHeuristic;
	const std::string sName;
};

const SHeuristic H_OTR(nmHeuristics::OTR, "OTR");
const SHeuristic H_ORT(nmHeuristics::ORT, "ORT");
const SHeuristic H_TRO(nmHeuristics::TRO, "TRO");
const SHeuristic H_TOR(nmHeuristics::TOR, "TOR");
const SHeuristic H_ROT(nmHeuristics::ROT, "ROT");
const SHeuristic H_RTO(nmHeuristics::RTO, "RTO");