#include "gtest/gtest.h"
#include "../MS_Solver/INC/Text.h"

TEST(Sentence, constructor_1) 
{
	Sentence s = "abc";
}
TEST(Sentence, constructor_2) 
{
	std::string str = "abc";
	Sentence s = str;
}
TEST(Sentence, constructor_3) 
{
	std::string str = "abc";
	Sentence s = str + "qwer";
}
TEST(Sentence, constructor_4)
{
	Sentence s = "qwer";
	Sentence s2 = s;
}
TEST(Sentence, constain_icase_1) 
{
	Sentence s = "cat, dog, bird, cow";
	EXPECT_TRUE(s.contain_icase("cat", "dog"));
}
TEST(Sentence, constain_icase_2)
{
	Sentence s = "cat, dog, bird, cow";
	EXPECT_FALSE(s.contain_icase("cat", "dog", "shrimp"));
}
TEST(Sentence, remove_after_1) {
	Sentence s = "abcdef";
	s.remove_after("d");

	Sentence ref = "abcd";
	EXPECT_EQ(s, ref);
}
TEST(Sentence, find_1) {
	Sentence s = "abcdef";
	const auto result = s.find_position("d");

	const auto ref = 3;
	EXPECT_EQ(result, ref);
}
TEST(Sentence, remove_from_here_1) {
	Sentence s = "abcdef";
	const auto pos = s.find_position("d");
	s.remove_from_here(pos);

	const Sentence ref = "abc";
	EXPECT_EQ(s, ref);
}
TEST(Sentence, remove_1) {
	Sentence s = "a1b2c1d2e1";
	s.remove({ '1','2' });

	Sentence ref = "abcde";
	EXPECT_EQ(s, ref);
}

TEST(Text, constructor_1) 
{
	Text txt = { "abc", "def", "ghi" };
}
TEST(Text, constructor_2)
{
	Sentence s = "qwer";
	Text txt = { s };
}
TEST(Text, merge_1) {
	Text txt = { "abc"};
	Text txt2 = { "def"};
	txt.merge(std::move(txt2));

	Text ref = { "abc", "def" };
	EXPECT_EQ(txt, ref);
}
TEST(Text, write_1) {
	Text txt;
	txt << "a \n\nb \n\nc ";
	txt.write("test.txt");

	Text txt2;
	txt2.read("test.txt");
	const auto result = txt2.size();

	const auto ref = 5;
	EXPECT_EQ(result, ref);
}
TEST(Text, add_empty_line_1)
{
	Text txt;
	txt << "test";
	txt << "test";
	txt << "test";
	txt.add_empty_lines(3);
	const auto result = txt.size();

	constexpr auto ref = 6;
	EXPECT_EQ(result, ref);
}

TEST(ms, contain_icase_1) 
{
	std::string str = "abc_qwer,wer__,,";
	std::string target_str = "abc";
	EXPECT_TRUE(ms::contains_icase(str, target_str));
}
TEST(ms, replace_1) {
	std::string str = "abc_qwer,wer__,,";
	ms::replace(str, ",", "_");

	std::string ref = "abc_qwer_wer____";
	EXPECT_EQ(str, ref);
}
TEST(ms, replace_2) {
	std::string str = "abc_qwer,wer__,,";
	ms::replace(str, "", "_");

	std::string ref = "abc_qwer,wer__,,";
	EXPECT_EQ(str, ref);
}
TEST(ms, replace_3) {
	std::string str = "abc_qwer,wer__,,";
	ms::replace(str, "wer", "");

	std::string ref = "abc_q,__,,";
	EXPECT_EQ(str, ref);
}
TEST(ms, replace_4) {
	std::string str = "abc_qwer,wer__,,";
	ms::replace(str, ',', '_');

	std::string ref = "abc_qwer_wer____";
	EXPECT_EQ(str, ref);
}
TEST(ms, replace_5) {
	std::string str = " target";
	ms::replace(str, " ", "");

	std::string ref = "target";
	EXPECT_EQ(str, ref);
}
TEST(ms, parse_1) {
	std::string str = "abc_qwer,wer__,,";
	const auto result = ms::parse(str, { ',', '_' });

	std::vector<std::string> ref = { "abc","qwer","wer" };
	EXPECT_EQ(result, ref);
}
TEST(ms, get_upper_case_1) {
	std::string str = "abc";
	const auto result = ms::get_upper_case(str);

	std::string ref = "ABC";
	EXPECT_EQ(result, ref);
}
TEST(ms, get_upper_case_2) {
	std::string str = "abc123";
	const auto result = ms::get_upper_case(str);

	std::string ref = "ABC123";
	EXPECT_EQ(result, ref);
}
TEST(ms, get_upper_case_3) {
	std::string str = "abc_123q";
	const auto result = ms::get_upper_case(str);

	std::string ref = "ABC_123Q";
	EXPECT_EQ(result, ref);
}
TEST(ms, get_upper_case_4) {
	const auto result = ms::get_upper_case("abc");
	std::string_view ref = "ABC";
	EXPECT_EQ(result, ref);
}
TEST(ms, find_icase_1) {
	std::string str = "abc_123q";
	const auto result = ms::find_icase(str, "BC_");

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
TEST(ms, find_icase_2) {
	std::string str = "abc_123q";
	const auto result = ms::find_icase(str, "C_12");

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}
TEST(ms, is_there_icase_1) {
	std::string str = "abc_123q";
	const auto result = ms::contains_icase(str, "C_12");

	const auto ref = true;
	EXPECT_EQ(result, ref);
}
TEST(ms, is_there_icase_2) {
	std::string str = "abc_123q";
	const auto result = ms::contains_icase(str, "3Q");

	const auto ref = true;
	EXPECT_EQ(result, ref);
}
TEST(ms, is_there_icase_3) {
	std::string str = "abc_123q";
	const auto result = ms::contains_icase(str, "a","q");

	const auto ref = true;
	EXPECT_EQ(result, ref);
}
TEST(ms, is_there_icase_4) {
	std::string str = "abc_123q";
	const auto result = ms::contains_icase(str, "ad", "q");

	const auto ref = false;
	EXPECT_EQ(result, ref);
}
TEST(ms, is_there_icase_5) {
	std::string str = "abc_123q";
	const auto result = ms::contains_icase(str, "ad", "qq");

	const auto ref = false;
	EXPECT_EQ(result, ref);
}
TEST(ms, rfind_nth_1) {
	std::string str = "abcaba";
	std::string target_str = "a";
	const auto result = ms::rfind_nth(str, target_str, 0);

	const auto ref = std::string::npos;
	EXPECT_EQ(result, ref);
}
TEST(ms, rfind_nth_2) {
	std::string str = "abcaba";
	std::string target_str = "a";
	const auto result = ms::rfind_nth(str, target_str, 1);

	const auto ref = 5;
	EXPECT_EQ(result, ref);
}
TEST(ms, rfind_nth_3) {
	std::string str = "abcaba";
	std::string target_str = "a";
	const auto result = ms::rfind_nth(str, target_str, 2);

	const auto ref = 3;
	EXPECT_EQ(result, ref);
}
TEST(ms, rfind_nth_4) {
	std::string str = "abcaba";
	std::string target_str = "a";
	const auto result = ms::rfind_nth(str, target_str, 3);

	const auto ref = 0;
	EXPECT_EQ(result, ref);
}
TEST(ms, rfind_nth_5) {
	std::string str = "abcaba";
	std::string target_str = "a";
	const auto result = ms::rfind_nth(str, target_str, 4);

	const auto ref = std::string::npos;
	EXPECT_EQ(result, ref);
}
TEST(ms, rfind_nth_6) {
	std::string str = "abcaba";
	std::string target_str = "b";
	const auto result = ms::rfind_nth(str, target_str, 3);

	const auto ref = std::string::npos;
	EXPECT_EQ(result, ref);
}
TEST(ms, remove_1) {
	std::string str = " target";
	ms::remove(str, " ");

	std::string ref = "target";
	EXPECT_EQ(str, ref);
}


//TEST(ms, rename) {
//	const std::string path = "D:/CODE/MS_Solver/MS_Solver/RSC/Grid/3D/";
//	const auto file_paths = ms::file_paths_in_path(path);
//
//	for (const auto& file_path : file_paths) {
//		const auto file_name = ms::remove(file_path, path);
//		const auto new_name = ms::replace(file_name, "Quad", "Hexa");
//		ms::rename(path, file_name, new_name);
//	}
//}

//#include "../MS_Solver/INC/Euclidean_Vector.h"
//TEST(ms, extract_file_name_1) {
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

