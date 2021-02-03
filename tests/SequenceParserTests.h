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

// New
class RawAdvMSATest : public ::testing::Test {
protected:
  void SetUp() override {
  }

};

TEST(TestRawAdvMSA, EqualEmpty) {
  IO::RawAdvMSA MSA1 = {0, 0, {}};
  IO::RawAdvMSA MSA2 = {0, 0, {}};
  EXPECT_EQ(MSA1, MSA2);
}

TEST(TestRawAdvMSA, EqualBasic) {
  IO::RawAdvMSA MSA1 = {2, 2, {}};
  MSA1.seqs["Taxa1"] = {{IO::StateFreq({'A', 1.0})}, {IO::StateFreq({'T', 1.0})}};
  MSA1.seqs["Taxa2"] = {{IO::StateFreq({'A', 0.8}), IO::StateFreq({'C', 0.2})},
			{IO::StateFreq({'T', 0.3}), IO::StateFreq({'G', 0.7})}};

  IO::RawAdvMSA MSA2 = {2, 2, {}};
  MSA2.seqs["Taxa1"] = {{IO::StateFreq({'A', 1.0})}, {IO::StateFreq({'T', 1.0})}};
  MSA2.seqs["Taxa2"] = {{IO::StateFreq({'C', 0.2}), IO::StateFreq({'A', 0.8})},
			{IO::StateFreq({'G', 0.7}), IO::StateFreq({'T', 0.3})}};


  EXPECT_EQ(MSA1, MSA2);
}

TEST(TestRawAdvMSA, NotEqual) {
  IO::RawAdvMSA MSA1 = {2, 2, {}};
  MSA1.seqs["Taxa1"] = {{IO::StateFreq({'A', 1.0})}, {IO::StateFreq({'T', 1.0})}};
  MSA1.seqs["Taxa2"] = {{IO::StateFreq({'G', 0.8}), IO::StateFreq({'C', 0.2})},
			{IO::StateFreq({'T', 0.3}), IO::StateFreq({'G', 0.7})}};

  IO::RawAdvMSA MSA2 = {2, 2, {}};
  MSA2.seqs["Taxa1"] = {{IO::StateFreq({'A', 1.0})}, {IO::StateFreq({'T', 1.0})}};
  MSA2.seqs["Taxa2"] = {{IO::StateFreq({'C', 0.2}), IO::StateFreq({'A', 0.8})},
			{IO::StateFreq({'G', 0.7}), IO::StateFreq({'T', 0.3})}};


  EXPECT_FALSE(MSA1 == MSA2);
}

TEST(TestRawAdvMSA, EmptyFile) {
  std::string data = "";

  try {
    IO::RawAdvMSA msa = IO::readRawAdvMSA(data, {"A", "T", "C", "G"});
    FAIL() << "Expected IO::ParseException"; 
  } catch(IO::ParseException const &err) {
    EXPECT_EQ(err.what(),std::string("empty file"));
  } catch(...) {
    FAIL() << "Expected IO::ParseException";
  }
}

TEST(TestRawAdvMSA, EmptySequence) {
  std::string data = ">TaxaA\n\n>Taxa1\nATATAT";

  try {
    IO::RawAdvMSA msa = IO::readRawAdvMSA(data, {"A", "T", "C", "G"});
    FAIL() << "Expected IO::ParseException"; 
  } catch(IO::ParseException const &err) {
    EXPECT_EQ(err.what(),std::string("expecting sequence name, possible missing \">\""));
  } catch(...) {
    FAIL() << "Expected IO::ParseException";
  }
}

TEST(TestRawAdvMSA, UnequalLength) {
  std::string data = ">TaxaA\nAT\n>Taxa1\nATATAT";

  try {
    IO::RawAdvMSA msa = IO::readRawAdvMSA(data, {"A", "T", "C", "G"});
    FAIL() << "Expected IO::ParseException"; 
  } catch(IO::ParseException const &err) {
    EXPECT_EQ(err.what(),std::string("sequences in fasta file are not equal length"));
  } catch(...) {
    FAIL() << "Expected IO::ParseException";
  }
}

TEST(TestRawAdvMSA, BasicParse) {
  std::string data = ">Taxa1\nATA";
  IO::RawAdvMSA MSA = {1, 3, {}};
  MSA.seqs["Taxa1"] = {{IO::StateFreq({'A', 1.0})}, {IO::StateFreq({'T', 1.0})}, {IO::StateFreq({'A', 1.0})}};

  EXPECT_EQ(IO::readRawAdvMSA(data, {"A", "T", "G", "C"}), MSA);
}

TEST(TestRawAdvMSA, InvalidSequence) {
  std::string data = ">Taxa1\nZTA";

  try {
    IO::RawAdvMSA msa = IO::readRawAdvMSA(data, {"A", "T", "C", "G"});
    FAIL() << "Expected IO::ParseException"; 
  } catch(IO::ParseException const &err) {
    EXPECT_EQ(err.what(),std::string("sequence for Taxa1 contains unrecognized state"));
  } catch(...) {
    FAIL() << "Expected IO::ParseException";
  }
}


