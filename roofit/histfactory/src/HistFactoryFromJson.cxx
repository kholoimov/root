#include "RooStats/HistFactory/HistFactoryFromJson.h"
#include <nlohmann/json.hpp>


using json = nlohmann::json;

HistFactoryFromJson::HistFactoryFromJson(const char* name, const char* title)
{
    measurement = RooStats::HistFactory::Measurement(name, title);
}

void HistFactoryFromJson::ReadAndConvertJSON(const std::string& filename)
{
    std::ifstream file(filename);
    file >> j;

    auto observations = j["observations"];
    // std::cout << observations << std::endl;

    for (auto& obs : observations)
    {
        std::string name = obs["name"];
        auto data = obs["data"];
        std::cout << name << " " << data << std::endl;
        // double uncertainty = obs["uncertainty"];
        // measurement.SetObservation(name, value, uncertainty);
    }
}