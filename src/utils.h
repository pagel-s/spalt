/**
 * @file utils.h
 * @brief Utility functions for system operations and string manipulation
 * @author Sebastian
 * @date 2024
 */

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef UTILS_H
#define UTILS_H

/**
 * @brief Check if a program is available in the system PATH
 *
 * Uses the 'which' command to check if the specified program is available
 * and executable in the current system PATH.
 *
 * @param program Name of the program to check (e.g., "msms", "xtb")
 * @return true if program is available, false otherwise
 */
bool isAvailable(const std::string& program);

/**
 * @brief Execute a system command
 *
 * Runs the specified command using the system shell and returns the success status.
 * The command output is captured and any errors are handled gracefully.
 *
 * @param command Command string to execute
 * @return true if command executed successfully (exit code 0), false otherwise
 */
bool run_command(const std::string& command);

/**
 * @brief Split a string into tokens based on whitespace
 *
 * Splits the input string into a vector of tokens using any whitespace character
 * (spaces, tabs, newlines) as delimiters. Empty tokens are ignored.
 *
 * @param str Input string to split
 * @return Vector of string tokens
 *
 * @note Uses std::istringstream for robust whitespace handling
 */
std::vector<std::string> split_on_whitespace(const std::string& str);

/**
 * @brief Safely remove a file or directory
 *
 * Attempts to remove the specified file or directory path. If removal fails,
 * the error is logged but no exception is thrown, making this safe for cleanup
 * operations where failure is not critical.
 *
 * @param p Filesystem path to remove
 *
 * @note Errors are logged to std::cerr but do not cause program termination
 */
void safe_remove(const std::filesystem::path& p);

/**
 * @brief Load all lines from a text file
 *
 * Reads all lines from the specified file and returns them as a vector of strings.
 * Empty lines are included in the result. This is useful for parsing structured
 * text files like coordinate files or data files.
 *
 * @param file Path to the file to load
 * @return Vector of strings, one per line in the file
 *
 * @throws std::runtime_error if file cannot be opened
 *
 * @note Lines are read using std::getline, so newline characters are stripped
 */
std::vector<std::string> load_file(std::string_view file);
#endif