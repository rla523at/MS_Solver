#include "gtest/gtest.h"
#include "../MS_Solver/INC/Text.h"

TEST(Text, merge_1) {
	Text txt = { "abc"};
	Text txt2 = { "def"};
	txt.merge(std::move(txt2));

	Text ref = { "abc", "def" };
	EXPECT_EQ(txt, ref);
}

GTEST_TEST(Text, write) {
	Text txt;
	txt << "a \n\nb \n\nc ";
	txt.write("test.txt");

	Text txt2;
	txt2.read("test.txt");
	const auto result = txt2.size();

	const auto ref = 5;
	EXPECT_EQ(result, ref);
}

GTEST_TEST(ms, replace_all_1) {
	std::string str = "abc_qwer,wer__,,";
	ms::be_replaced(str, ",", "_");

	std::string ref = "abc_qwer_wer____";
	EXPECT_EQ(str, ref);
}
GTEST_TEST(ms, replace_all_2) {
	std::string str = "abc_qwer,wer__,,";
	ms::be_replaced(str, "", "_");

	std::string ref = "abc_qwer,wer__,,";
	EXPECT_EQ(str, ref);
}
GTEST_TEST(ms, replace_all_3) {
	std::string str = "abc_qwer,wer__,,";
	ms::be_replaced(str, "wer", "");

	std::string ref = "abc_q,__,,";
	EXPECT_EQ(str, ref);
}
GTEST_TEST(ms, replace_all_4) {
	std::string str = "abc_qwer,wer__,,";
	ms::be_replaced(str, ',', '_');

	std::string ref = "abc_qwer_wer____";
	EXPECT_EQ(str, ref);
}

GTEST_TEST(ms, parse_1) {
	std::string str = "abc_qwer,wer__,,";
	const auto result = ms::parse(str, { ',', '_' });

	std::vector<std::string> ref = { "abc","qwer","wer" };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(ms, upper_case_1) {
	std::string str = "abc";
	const auto result = ms::upper_case(str);

	std::string ref = "ABC";
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, upper_case_2) {
	std::string str = "abc123";
	const auto result = ms::upper_case(str);

	std::string ref = "ABC123";
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, upper_case_3) {
	std::string str = "abc_123q";
	const auto result = ms::upper_case(str);

	std::string ref = "ABC_123Q";
	EXPECT_EQ(result, ref);
}

GTEST_TEST(ms, find_icase_1) {
	std::string str = "abc_123q";
	const auto result = ms::find_icase(str, "BC_");

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, find_icase_2) {
	std::string str = "abc_123q";
	const auto result = ms::find_icase(str, "C_12");

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}

GTEST_TEST(ms, is_there_icase_1) {
	std::string str = "abc_123q";
	const auto result = ms::contains_icase(str, "C_12");

	const auto ref = true;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, is_there_icase_2) {
	std::string str = "abc_123q";
	const auto result = ms::contains_icase(str, "3Q");

	const auto ref = true;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, is_there_icase_3) {
	std::string str = "abc_123q";
	const auto result = ms::contains_icase(str, "a","q");

	const auto ref = true;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, is_there_icase_4) {
	std::string str = "abc_123q";
	const auto result = ms::contains_icase(str, "ad", "q");

	const auto ref = false;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, is_there_icase_5) {
	std::string str = "abc_123q";
	const auto result = ms::contains_icase(str, "ad", "qq");

	const auto ref = false;
	EXPECT_EQ(result, ref);
}


GTEST_TEST(ms, rfind_nth_1) {
	std::string str = "abcaba";
	std::string target_str = "a";
	const auto result = ms::rfind_nth(str, target_str, 0);

	const auto ref = std::string::npos;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, rfind_nth_2) {
	std::string str = "abcaba";
	std::string target_str = "a";
	const auto result = ms::rfind_nth(str, target_str, 1);

	const auto ref = 5;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, rfind_nth_3) {
	std::string str = "abcaba";
	std::string target_str = "a";
	const auto result = ms::rfind_nth(str, target_str, 2);

	const auto ref = 3;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, rfind_nth_4) {
	std::string str = "abcaba";
	std::string target_str = "a";
	const auto result = ms::rfind_nth(str, target_str, 3);

	const auto ref = 0;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, rfind_nth_5) {
	std::string str = "abcaba";
	std::string target_str = "a";
	const auto result = ms::rfind_nth(str, target_str, 4);

	const auto ref = std::string::npos;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(ms, rfind_nth_6) {
	std::string str = "abcaba";
	std::string target_str = "b";
	const auto result = ms::rfind_nth(str, target_str, 3);

	const auto ref = std::string::npos;
	EXPECT_EQ(result, ref);
}

TEST(ms, be_removed_1) {
	std::string str = " target";
	ms::be_removed(str, " ");

	std::string ref = "target";
	EXPECT_EQ(str, ref);
}

TEST(ms, be_replaced_1) {
	std::string str = " target";
	ms::be_replaced(str, " ", "");

	std::string ref = "target";
	EXPECT_EQ(str, ref);
}

//TEST(ms, rename) {
//	const std::string path = "D:/CODE/MS_Solver/MS_Solver/RSC/Grid/3D/";
//	const auto file_paths = ms::file_paths_in_path(path);
//
//	for (const auto& file_path : file_paths) {
//		const auto file_name = ms::remove(file_path, path);
//		const auto new_name = ms::replace_all(file_name, "Quad", "Hexa");
//		ms::rename(path, file_name, new_name);
//	}
//}

//#include "../MS_Solver/INC/Euclidean_Vector.h"
//GTEST_TEST(ms, extract_file_name_1) {
//	const std::string path = "C:/Users/KimMinSeok/source/repos/MS_Test/MS_Test/RSC/Quadrature/Standard/Quadrilateral/";
//	const auto path_txt = ms::extract_file_path_text(path);
//
//	for (const auto& file_path : path_txt) {
//		Text quadrature_txt;
//		quadrature_txt.read_line_by_line(file_path);
//
//		std::string node_str = "{ ";
//		std::string weight_str = "{ ";
//		for (const auto& sentence : quadrature_txt) {
//			const char denominator = ' ';
//			const auto parsed_sentence = ms::parse(sentence, denominator);
//			node_str += "{ " + parsed_sentence[0] + ", " + parsed_sentence[1] + " }, ";
//			weight_str += parsed_sentence[2] + ", ";
//		}
//		node_str.pop_back();
//		node_str.pop_back();
//		weight_str.pop_back();
//		weight_str.pop_back();
//
//		node_str += " }";
//		weight_str += " }";
//
//		const auto file_name = ms::erase(file_path, path);
//		const auto parsed_strs = ms::parse(file_name, '_');
//		const auto& order_str = parsed_strs[0];
//
//		auto quadrature_str = "{ " + node_str + ", " + weight_str + " }";
//		Text modified_txt = { quadrature_str };
//		modified_txt.write(path + "modified/" + order_str + ".txt");
//	}
//
//	//std::cout << name_txt;
//}

