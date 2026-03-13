#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// Forward declaration of the function we want to test from cli.cpp
std::vector<std::pair<std::string, std::string>> loadSmilesFromFile(const std::string& file_path);

class TSVParsingTest : public ::testing::Test {
  protected:
    void SetUp() override {
        test_file_path = "test_input.txt";
    }

    void TearDown() override {
        if (std::filesystem::exists(test_file_path)) {
            std::filesystem::remove(test_file_path);
        }
    }

    void CreateTestFile(const std::string& content) {
        std::ofstream file(test_file_path);
        file << content;
        file.close();
    }

    std::string test_file_path;
};

TEST_F(TSVParsingTest, ParseSingleColumn) {
    CreateTestFile("CCO\nc1ccccc1\n");
    auto result = loadSmilesFromFile(test_file_path);
    ASSERT_EQ(result.size(), 2);
    EXPECT_EQ(result[0].first, "CCO");
    EXPECT_EQ(result[0].second, "");
    EXPECT_EQ(result[1].first, "c1ccccc1");
    EXPECT_EQ(result[1].second, "");
}

TEST_F(TSVParsingTest, ParseTabSeparated) {
    CreateTestFile("CCO\tEthanol\nc1ccccc1\tBenzene\n");
    auto result = loadSmilesFromFile(test_file_path);
    ASSERT_EQ(result.size(), 2);
    EXPECT_EQ(result[0].first, "CCO");
    EXPECT_EQ(result[0].second, "Ethanol");
    EXPECT_EQ(result[1].first, "c1ccccc1");
    EXPECT_EQ(result[1].second, "Benzene");
}

TEST_F(TSVParsingTest, ParseSpaceSeparated) {
    CreateTestFile("CCO Ethanol\nc1ccccc1 Benzene\n");
    auto result = loadSmilesFromFile(test_file_path);
    ASSERT_EQ(result.size(), 2);
    EXPECT_EQ(result[0].first, "CCO");
    EXPECT_EQ(result[0].second, "Ethanol");
    EXPECT_EQ(result[1].first, "c1ccccc1");
    EXPECT_EQ(result[1].second, "Benzene");
}

TEST_F(TSVParsingTest, ParseWithWhitespace) {
    CreateTestFile("  CCO \t Ethanol  \n\tc1ccccc1   Benzene  \n");
    auto result = loadSmilesFromFile(test_file_path);
    ASSERT_EQ(result.size(), 2);
    EXPECT_EQ(result[0].first, "CCO");
    EXPECT_EQ(result[0].second, "Ethanol");
    EXPECT_EQ(result[1].first, "c1ccccc1");
    EXPECT_EQ(result[1].second, "Benzene");
}

TEST_F(TSVParsingTest, ParseWithComments) {
    CreateTestFile("# This is a comment\nCCO\tEthanol\n# Another comment\nc1ccccc1\tBenzene\n");
    auto result = loadSmilesFromFile(test_file_path);
    ASSERT_EQ(result.size(), 2);
    EXPECT_EQ(result[0].first, "CCO");
    EXPECT_EQ(result[0].second, "Ethanol");
    EXPECT_EQ(result[1].first, "c1ccccc1");
    EXPECT_EQ(result[1].second, "Benzene");
}
