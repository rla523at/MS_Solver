

//TEST(VECTORFUNCTION, OPERATOR_CALL1) {
//	auto p1 = (X ^ 2) + 3 * (X ^ 2) * Y + (Y ^ 3) + (Z ^ 2) - 6;
//	auto p2 = X + Y + Z - 3;
//	auto p3 = (Y ^ 2) * Z + X * Z - 2;
//	VectorFunction<auto> f = { p1,p2,p3 };
//	MathVector node = { 1,1,1 };
//	const auto result = f(node);
//
//	MathVector ref = { 0,0,0 };
//	EXPECT_EQ(result, ref);
//}
//TEST(VECTORFUNCTION, DIFFERENTIATE1) {
//	auto p1 = X + 1;
//	auto p2 = X + Y + Z - 3;
//	auto p3 = Y + Z;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate(variable_index);
//
//	VectorFunction<auto> ref = { 1,1,0 };
//	EXPECT_EQ(f, ref);
//}
//TEST(VECTORFUNCTION, DIFFERENTIATE2) {
//	auto p1 = X * Y + 1;
//	auto p2 = X + Y * Z + Z - 3;
//	auto p3 = Y + Z;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate(variable_index);
//
//	VectorFunction<auto> ref = { Y,1,0 };
//	EXPECT_EQ(f, ref);
//}
//TEST(VECTORFUNCTION, DIFFERENTIATE3) {
//	auto p1 = X * Y + 1;
//	auto p2 = X + Y * Z + Z - 3;
//	auto p3 = Y + Z;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 1;
//	f.differentiate(variable_index);
//
//	VectorFunction<auto> ref = { X,Z,1 };
//	EXPECT_EQ(f, ref);
//}
//TEST(VECTORFUNCTION, DIFFERENTIATE4) {
//	auto p1 = 1.5 * X + 0.5 * Y + 3;
//	auto p2 = Y + 3;
//	auto p3 = 0;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate(variable_index);
//
//	VectorFunction<auto> ref = { 1.5,0,0 };
//	EXPECT_EQ(f, ref);
//}
//TEST(VECTORFUNCTION, DIFFERENTIATE5) {
//	auto p1 = 1.5 * X + 0.5 * Y + 3;
//	auto p2 = Y + 3;
//	auto p3 = 0;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 1;
//	f.differentiate(variable_index);
//
//	VectorFunction<auto> ref = { 0.5,1,0 };
//	EXPECT_EQ(f, ref);
//}
//TEST(VECTORFUNCTION, CROSS_PRODUCT1) {
//	VectorFunction<auto> vf1 = { 1.5,0,0 };
//	VectorFunction<auto> vf2 = { 0.5,1,0 };
//	const auto result = vf1.cross_product(vf2);
//
//	VectorFunction<auto> ref = { 0,0,1.5 };
//	EXPECT_EQ(result, ref);
//}
//
//
//TEST(VECTORFUNCTION, L2_NORM1) {
//	VectorFunction<auto> vf = { 0,0,1.5 };
//	const auto result = vf.L2_norm();
//
//	auto ref = 1.5;
//	EXPECT_EQ(result, ref);
//}
