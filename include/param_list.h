#ifndef PARAM_LIST
#define PARAM_LIST

#include <string>
#include <unordered_map>
#include <variant>

typedef std::unordered_map<std::string, std::variant<bool, int, double>> param_list;

#endif
