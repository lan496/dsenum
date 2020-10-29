#include <iostream>
#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>
#include <permutation.hpp>
#include <type.hpp>
#include <spec/isomorphism.hpp>
using namespace permutation;
using namespace tdzdd;

void test1() {
    auto perm = Permutation(std::vector<Element>{2, 1, 0});
    PermutationFrontierManager pfm(perm);
    pfm.dump(std::cerr);

    IsomorphismElimination spec(pfm);

    tdzdd::MessageHandler mh;
    mh.begin("begin");

    DdStructure<2> dd(spec);
    dd.zddReduce();

    mh.end();

    auto actual = dd.zddCardinality();
    std::cerr << "# of solutions: " << actual << std::endl;

    std::ofstream output("debug.dot");
    dd.dumpDot(output);

    assert(actual == "6");
}

int main() {
    tdzdd::MessageHandler::showMessages(true);
    test1();

    return 0;
}
