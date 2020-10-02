#include <iostream>
#include "../src/IO/SequencesParser.h"

class RawMSATest : public ::testing::Test {
 protected:
  void SetUp() override {
    base_msa.n = 2;
    base_msa.cols = 8;
    base_msa.seqs["taxaA"] = "ATCGATCG";
    base_msa.seqs["taxaB"] = "GCTAGCTA";

    gap_msa.n = 2;
    gap_msa.cols = 8;
    gap_msa.seqs["taxaA"] = "-TCG-TCG";
    gap_msa.seqs["taxaB"] = "GCT-GCT-";
  }

  IO::RawMSA base_msa = {};
  IO::RawMSA gap_msa = {};
};

TEST_F(RawMSATest, BasicFunctionality) {
  std::list<std::string> remove_list = {"A"};
  IO::convertToGaps(base_msa, remove_list);

  EXPECT_EQ(base_msa, gap_msa);
  EXPECT_EQ(IO::getRawMSANames(base_msa), std::list<std::string>({"taxaA", "taxaB"}));
}

