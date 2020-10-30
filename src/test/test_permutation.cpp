#include <iostream>
#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>
#include <permutation.hpp>
#include <type.hpp>
#include <spec/isomorphism.hpp>
using namespace permutation;
using namespace tdzdd;

void test1(bool dump_dot) {
    std::vector<Element> sigma = {2, 1, 0};
    auto perm = Permutation(sigma);
    PermutationFrontierManager pfm(perm);
    // pfm.dump(std::cerr);

    isomorphism::IsomorphismElimination spec(pfm);

    tdzdd::MessageHandler mh;
    mh.begin("begin");

    DdStructure<2> dd(spec);
    dd.zddReduce();

    mh.end();

    auto actual = dd.zddCardinality();
    std::cout << "# of solutions: " << actual << std::endl;

    if (dump_dot) {
        std::ofstream output("debug.dot");
        dd.dumpDot(output);
    }

    assert(actual == "6");
}

void test2(bool dump_dot) {
    auto fname_dot = "debug2.dot";
    auto perm = Permutation(std::vector<Element>{1, 2, 0});
    std::string cardinarlity_expect = "5";

    PermutationFrontierManager pfm(perm);
    pfm.dump(std::cerr);

    isomorphism::IsomorphismElimination spec(pfm);

    tdzdd::MessageHandler mh;
    mh.begin("begin");

    DdStructure<2> dd(spec);
    // dd.zddReduce();

    mh.end();

    auto actual = dd.zddCardinality();
    std::cout << "# of solutions: " << actual << std::endl;

    if (dump_dot) {
        std::ofstream output(fname_dot);
        dd.dumpDot(output);
    }

    assert(actual == cardinarlity_expect);
}

int main() {
    tdzdd::MessageHandler::showMessages(true);
    test1(false);
    test2(true);

    return 0;
}
