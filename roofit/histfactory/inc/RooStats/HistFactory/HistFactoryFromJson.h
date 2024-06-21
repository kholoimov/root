#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include "Measurement.h"
#include "Channel.h"
#include "Sample.h"

using json = nlohmann::json;

#pragma once

class HistFactoryFromJson {
public:
    HistFactoryFromJson(const char* name, const char* title);

    ~HistFactoryFromJson() {}

    void ReadAndConvertJSON(const std::string& filename);

    RooStats::HistFactory::Measurement getMeasurement() 
    {   
        return measurement;
    }

private:
    json j;
    RooStats::HistFactory::Measurement measurement;
};