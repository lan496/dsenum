#include <iostream>
#include <string>
#include <unordered_set>
#include <cassert>

#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>

#include <type.hpp>
#include <spec/universe.hpp>

using namespace pyzdd;
using namespace pyzdd::universe;

void check_enumerated(int n, const std::string& expect) {
    Universe spec(n);
    tdzdd::DdStructure<2> dd(spec);
#ifndef _DEBUG
    dd.zddReduce();
#endif

    auto actual = dd.zddCardinality();
#ifdef _DEBUG
    std::cerr << "# of solutions: " << actual << std::endl;
    std::ofstream output("debug.dot");
    dd.dumpDot(output);
#endif

    assert(actual == expect);
}

void test1() {
    int n = 4;
    std::string cardinality_expect = "16";
    check_enumerated(n, cardinality_expect);
}

int main() {
    tdzdd::MessageHandler::showMessages(true);
    check_enumerated(4, "16");
    check_enumerated(10, "1024");
    return 0;
}
