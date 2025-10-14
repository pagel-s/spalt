#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "utils.h"

bool isAvailable(const std::string& program) {
    if (const char* path = std::getenv("PATH")) {
        std::string p(path);
        size_t pos = 0;
        while ((pos = p.find(':')) != std::string::npos) {
            auto dir = p.substr(0, pos);
            if (std::filesystem::exists(dir + "/" + program)) {
                return true;
            }
            p.erase(0, pos + 1);
        }
    }
    return false;
}

bool run_command(const std::string& command) {
    int rc = std::system(command.c_str());
    if (rc == -1)
        return false;
    return WIFEXITED(rc) && WEXITSTATUS(rc) == 0;
}

std::vector<std::string> split_on_whitespace(const std::string& str) {
    std::vector<std::string> result;
    std::istringstream iss(str);
    std::string token;
    while (iss >> token) {
        result.push_back(token);
    }
    return result;
}

void safe_remove(const std::filesystem::path& p) {
    std::error_code ec;
    bool removed = std::filesystem::remove(p, ec);

    if (ec) {
#ifdef DEBUG
        std::cerr << "Error deleting " << p << ": " << ec.message() << "\n";
#endif
    } else if (!removed) {
#ifdef DEBUG
        std::cerr << "File not found: " << p << "\n";
#endif
    } else {
#ifdef DEBUG
        std::cout << "Deleted: " << p << "\n";
#endif
    }
}

std::vector<std::string> load_file(std::string_view file) {
    std::ifstream ifs{std::string(file)};
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line)) {
        lines.push_back(line);
    }
    return lines;
}