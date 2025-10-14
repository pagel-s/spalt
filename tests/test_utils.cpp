#include <gtest/gtest.h>
#include <fstream>
#include "utils.h"

TEST(Utils, SplitOnWhitespace) {
    std::string str = "1 2 3 4 5";
    std::vector<std::string> result = split_on_whitespace(str);
    EXPECT_EQ(result.size(), 5);
    EXPECT_EQ(result[0], "1");
    EXPECT_EQ(result[1], "2");

    std::string str2 = "27.842    0.076   34.073 2.0000 1 A_1_x_PRO_CA_x";
    std::vector<std::string> result2 = split_on_whitespace(str2);
    EXPECT_EQ(result2.size(), 6);
    EXPECT_EQ(result2[0], "27.842");
    EXPECT_EQ(result2[1], "0.076");
    EXPECT_EQ(result2[2], "34.073");
    EXPECT_EQ(result2[3], "2.0000");
    EXPECT_EQ(result2[4], "1");
    EXPECT_EQ(result2[5], "A_1_x_PRO_CA_x");
}

TEST(Utils, LoadFile) {
    std::ofstream file("test.esp");
    file << "1\n2\n3\n4\n5\n";
    file.close();

    std::vector<std::string> result = load_file("test.esp");
    EXPECT_EQ(result.size(), 5);
    EXPECT_EQ(result[0], "1");
    EXPECT_EQ(result[1], "2");
    EXPECT_EQ(result[2], "3");
    EXPECT_EQ(result[3], "4");
    EXPECT_EQ(result[4], "5");

    safe_remove("test.esp");
}