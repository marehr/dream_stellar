
#include <gtest/gtest.h>

#include <stellar/utils/fraction.hpp>

TEST(fraction, construct)
{
    {
        stellar::utils::fraction frac{0, 1};
        EXPECT_EQ(frac.numerator(), 0);
        EXPECT_EQ(frac.denominator(), 1);
    }

    {
        stellar::utils::fraction frac{15, 12};
        EXPECT_EQ(frac.numerator(), 15);
        EXPECT_EQ(frac.denominator(), 12);
    }

    {
        stellar::utils::fraction frac{12, 15};
        EXPECT_EQ(frac.numerator(), 12);
        EXPECT_EQ(frac.denominator(), 15);
    }

    {
        stellar::utils::fraction frac{-12, 15};
        EXPECT_EQ(frac.numerator(), -12);
        EXPECT_EQ(frac.denominator(), 15);
    }

    {
        stellar::utils::fraction frac{12, -15};
        EXPECT_EQ(frac.numerator(), -12);
        EXPECT_EQ(frac.denominator(), 15);
    }

    {
        stellar::utils::fraction frac{-12, -15};
        EXPECT_EQ(frac.numerator(), 12);
        EXPECT_EQ(frac.denominator(), 15);
    }
}

TEST(fraction, inverse)
{
    {
        EXPECT_THROW(inverse(stellar::utils::fraction{0, 1}), std::domain_error);
    }

    {
        stellar::utils::fraction frac = inverse(stellar::utils::fraction{15, 12});
        EXPECT_EQ(frac.numerator(), 12);
        EXPECT_EQ(frac.denominator(), 15);
    }

    {
        stellar::utils::fraction frac = inverse(stellar::utils::fraction{12, 15});
        EXPECT_EQ(frac.numerator(), 15);
        EXPECT_EQ(frac.denominator(), 12);
    }

    {
        stellar::utils::fraction frac = inverse(stellar::utils::fraction{-12, 15});
        EXPECT_EQ(frac.numerator(), -15);
        EXPECT_EQ(frac.denominator(), 12);
    }

    {
        stellar::utils::fraction frac = inverse(stellar::utils::fraction{12, -15});
        EXPECT_EQ(frac.numerator(), -15);
        EXPECT_EQ(frac.denominator(), 12);
    }

    {
        stellar::utils::fraction frac = inverse(stellar::utils::fraction{-12, -15});
        EXPECT_EQ(frac.numerator(), 15);
        EXPECT_EQ(frac.denominator(), 12);
    }
}

TEST(fraction, is_proper)
{
    EXPECT_TRUE((stellar::utils::fraction{0, 1}.is_proper()));
    EXPECT_TRUE((stellar::utils::fraction{12, 15}.is_proper()));
    EXPECT_TRUE((stellar::utils::fraction{-12, 15}.is_proper()));
    EXPECT_TRUE((stellar::utils::fraction{12, -15}.is_proper()));
    EXPECT_TRUE((stellar::utils::fraction{-12, -15}.is_proper()));

    EXPECT_FALSE((stellar::utils::fraction{15, 12}.is_proper()));
    EXPECT_FALSE((stellar::utils::fraction{-15, 12}.is_proper()));
    EXPECT_FALSE((stellar::utils::fraction{15, -12}.is_proper()));
    EXPECT_FALSE((stellar::utils::fraction{-15, -12}.is_proper()));
}

TEST(fraction, double)
{
    EXPECT_FLOAT_EQ((double)(stellar::utils::fraction{0, 1}), 0.0 / 1);
    EXPECT_FLOAT_EQ((double)(stellar::utils::fraction{12, 15}), 12.0 / 15);
    EXPECT_FLOAT_EQ((double)(stellar::utils::fraction{-12, 15}), -12.0 / 15);
    EXPECT_FLOAT_EQ((double)(stellar::utils::fraction{12, -15}), 12.0 / -15);
    EXPECT_FLOAT_EQ((double)(stellar::utils::fraction{-12, -15}), -12.0 / -15);

    EXPECT_FLOAT_EQ((double)(stellar::utils::fraction{15, 12}), 15.0 / 12);
    EXPECT_FLOAT_EQ((double)(stellar::utils::fraction{-15, 12}), -15.0 / 12);
    EXPECT_FLOAT_EQ((double)(stellar::utils::fraction{15, -12}), 15.0 / -12);
    EXPECT_FLOAT_EQ((double)(stellar::utils::fraction{-15, -12}), -15.0 / -12);
}

TEST(fraction, minus)
{
    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{0, 1};
        stellar::utils::fraction frac2 = stellar::utils::fraction{0, 10};
        stellar::utils::fraction frac = frac1 - frac2;
        EXPECT_EQ(frac.numerator(), 0);
        EXPECT_EQ(frac.denominator(), 10);
    }

    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{15, 12};
        stellar::utils::fraction frac2 = stellar::utils::fraction{15, 12};
        stellar::utils::fraction frac = frac1 - frac2;
        EXPECT_EQ(frac.numerator(), 0);
        EXPECT_EQ(frac.denominator(), 144);
    }

    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{12, 15};
        stellar::utils::fraction frac2 = stellar::utils::fraction{4, 5};
        stellar::utils::fraction frac = frac1 - frac2;
        EXPECT_EQ(frac.numerator(), 0);
        EXPECT_EQ(frac.denominator(), 15 * 5);
    }

    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{-12, 15};
        stellar::utils::fraction frac2 = stellar::utils::fraction{4, 5};
        stellar::utils::fraction frac = frac1 - frac2;
        EXPECT_EQ(frac.numerator(), -12 * 5 - 4 * 15);
        EXPECT_EQ(frac.denominator(), 15 * 5);
    }

    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{12, -15};
        stellar::utils::fraction frac2 = stellar::utils::fraction{4, 5};
        stellar::utils::fraction frac = frac1 - frac2;
        EXPECT_EQ(frac.numerator(), -12 * 5 - 4 * 15);
        EXPECT_EQ(frac.denominator(), 15 * 5);
    }

    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{-12, -15};
        stellar::utils::fraction frac2 = stellar::utils::fraction{4, 5};
        stellar::utils::fraction frac = frac1 - frac2;
        EXPECT_EQ(frac.numerator(), 12 * 5 - 4 * 15);
        EXPECT_EQ(frac.denominator(), 15 * 5);
    }
}

TEST(fraction, division)
{
    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{0, 1};
        stellar::utils::fraction frac2 = stellar::utils::fraction{0, 10};
        {
            EXPECT_THROW(frac1 / frac2, std::domain_error);
        }
    }

    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{15, 12};
        stellar::utils::fraction frac2 = stellar::utils::fraction{15, 12};
        stellar::utils::fraction frac = frac1 / frac2;
        EXPECT_EQ(frac.numerator(), 15 * 12);
        EXPECT_EQ(frac.denominator(), 12 * 15);
    }

    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{12, 15};
        stellar::utils::fraction frac2 = stellar::utils::fraction{4, 5};
        stellar::utils::fraction frac = frac1 / frac2;
        EXPECT_EQ(frac.numerator(), 12 * 5);
        EXPECT_EQ(frac.denominator(), 15 * 4);
    }

    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{-12, 15};
        stellar::utils::fraction frac2 = stellar::utils::fraction{4, 5};
        stellar::utils::fraction frac = frac1 / frac2;
        EXPECT_EQ(frac.numerator(), -12 * 5);
        EXPECT_EQ(frac.denominator(), 15 * 4);
    }

    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{12, -15};
        stellar::utils::fraction frac2 = stellar::utils::fraction{4, 5};
        stellar::utils::fraction frac = frac1 / frac2;
        EXPECT_EQ(frac.numerator(), - (12 * 5));
        EXPECT_EQ(frac.denominator(), std::abs(-15) * 4);
    }

    {
        stellar::utils::fraction frac1 = stellar::utils::fraction{-12, -15};
        stellar::utils::fraction frac2 = stellar::utils::fraction{4, 5};
        stellar::utils::fraction frac = frac1 / frac2;
        EXPECT_EQ(frac.numerator(), 12 * 5);
        EXPECT_EQ(frac.denominator(), 15 * 4);
    }
}

TEST(fraction, floor)
{
    EXPECT_EQ(floor(stellar::utils::fraction{0, 1}), 0); // ⌊0⌋ = 0
    EXPECT_EQ(floor(stellar::utils::fraction{1, 1}), 1); // ⌊1⌋ = 1
    EXPECT_EQ(floor(stellar::utils::fraction{2, 2}), 1); // ⌊1⌋ = -1
    EXPECT_EQ(floor(stellar::utils::fraction{-1, 1}), -1); // ⌊-1⌋ = -1
    EXPECT_EQ(floor(stellar::utils::fraction{-2, 2}), -1); // ⌊-1⌋ = -1
    EXPECT_EQ(floor(stellar::utils::fraction{1, 2}), 0); // ⌊0.5⌋ = 0
    EXPECT_EQ(floor(stellar::utils::fraction{-1, 2}), -1); // ⌊-0.5⌋ = -1
    EXPECT_EQ(floor(stellar::utils::fraction{5, 2}), 2); // ⌊2.5⌋ = 2
    EXPECT_EQ(floor(stellar::utils::fraction{-5, 2}), -3); // ⌊-2.5⌋ = -3
}

TEST(fraction, ceil)
{
    EXPECT_EQ(ceil(stellar::utils::fraction{0, 1}), 0); // ⌈0⌉ = 0
    EXPECT_EQ(ceil(stellar::utils::fraction{1, 1}), 1); // ⌈1⌉ = 1
    EXPECT_EQ(ceil(stellar::utils::fraction{2, 2}), 1); // ⌈1⌉ = 1
    EXPECT_EQ(ceil(stellar::utils::fraction{-1, 1}), -1); // ⌈-1⌉ = -1
    EXPECT_EQ(ceil(stellar::utils::fraction{-2, 2}), -1); // ⌈-1⌉ = -1
    EXPECT_EQ(ceil(stellar::utils::fraction{1, 2}), 1); // ⌈0.5⌉ = 1
    EXPECT_EQ(ceil(stellar::utils::fraction{-1, 2}), 0); // ⌈-0.5⌉ = 0
    EXPECT_EQ(ceil(stellar::utils::fraction{5, 2}), 3); // ⌈2.5⌉ = 3
    EXPECT_EQ(ceil(stellar::utils::fraction{-5, 2}), -2); // ⌈-2.5⌉ = -2
}

TEST(fraction, limit_denominator)
{
    {
        stellar::utils::fraction frac = stellar::utils::fraction{0, 1}.limit_denominator();
        EXPECT_EQ(frac.numerator(), 0);
        EXPECT_EQ(frac.denominator(), 1);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction{46116860184273879ll, 4611686018427387904ull}.limit_denominator();
        EXPECT_EQ(frac.numerator(), 1);
        EXPECT_EQ(frac.denominator(), 100);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction{50, 100}.limit_denominator(10);
        EXPECT_EQ(frac.numerator(), 1);
        EXPECT_EQ(frac.denominator(), 2);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction{100, 50}.limit_denominator(10);
        EXPECT_EQ(frac.numerator(), 2);
        EXPECT_EQ(frac.denominator(), 1);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction{2180466202566ll, 5451165506414999ull}.limit_denominator();
        EXPECT_EQ(frac.numerator(), 1);
        EXPECT_EQ(frac.denominator(), 2500);
    }

    // negative values
    {
        stellar::utils::fraction frac = stellar::utils::fraction{-2180466202566ll, 5451165506414999ull}.limit_denominator();
        EXPECT_EQ(frac.numerator(), -1);
        EXPECT_EQ(frac.denominator(), 2500);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction{-46116860184273879ll, 4611686018427387904ull}.limit_denominator();
        EXPECT_EQ(frac.numerator(), -1);
        EXPECT_EQ(frac.denominator(), 100);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction{-50, 100}.limit_denominator(10);
        EXPECT_EQ(frac.numerator(), -1);
        EXPECT_EQ(frac.denominator(), 2);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction{-100, 50}.limit_denominator(10);
        EXPECT_EQ(frac.numerator(), -2);
        EXPECT_EQ(frac.denominator(), 1);
    }
}

TEST(fraction, from_string)
{
    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0");
        EXPECT_EQ(frac.numerator(), 0);
        EXPECT_EQ(frac.denominator(), 1);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0.01");
        EXPECT_EQ(frac.numerator(), 1);
        EXPECT_EQ(frac.denominator(), 100);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0.2");
        EXPECT_EQ(frac.numerator(), 2);
        EXPECT_EQ(frac.denominator(), 10);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0.3");
        EXPECT_EQ(frac.numerator(), 3);
        EXPECT_EQ(frac.denominator(), 10);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0.4");
        EXPECT_EQ(frac.numerator(), 4);
        EXPECT_EQ(frac.denominator(), 10);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0.5");
        EXPECT_EQ(frac.numerator(), 5);
        EXPECT_EQ(frac.denominator(), 10);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0.6");
        EXPECT_EQ(frac.numerator(), 6);
        EXPECT_EQ(frac.denominator(), 10);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0.7");
        EXPECT_EQ(frac.numerator(), 7);
        EXPECT_EQ(frac.denominator(), 10);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0.8");
        EXPECT_EQ(frac.numerator(), 8);
        EXPECT_EQ(frac.denominator(), 10);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0.9");
        EXPECT_EQ(frac.numerator(), 9);
        EXPECT_EQ(frac.denominator(), 10);
    }

    {
        // note: uses double conversion, 0.0004
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("4e-4").limit_denominator();
        EXPECT_EQ(frac.numerator(), 1);
        EXPECT_EQ(frac.denominator(), 2500);
    }

    {
        // note: uses double conversion, -0.0004
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("-4e-4").limit_denominator();
        EXPECT_EQ(frac.numerator(), -1);
        EXPECT_EQ(frac.denominator(), 2500);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("0.000976562");
        EXPECT_EQ(frac.numerator(), 976'562);
        EXPECT_EQ(frac.denominator(), 1'000'000'000);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("1024");
        EXPECT_EQ(frac.numerator(), 1024);
        EXPECT_EQ(frac.denominator(), 1);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("1024.25");
        EXPECT_EQ(frac.numerator(), 102425);
        EXPECT_EQ(frac.denominator(), 100);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("999.999""999""999""99");
        EXPECT_EQ(frac.numerator(), 99'999'999'999'999);
        EXPECT_EQ(frac.denominator(), 100'000'000'000);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("999.999""999""999""999");
        EXPECT_EQ(frac.numerator(), 999'999'999'999'999);
        EXPECT_EQ(frac.denominator(), 1'000'000'000'000);
    }

    {
        stellar::utils::fraction frac = stellar::utils::fraction::from_string("1.8");
        EXPECT_EQ(frac.numerator(), 18);
        EXPECT_EQ(frac.denominator(), 10);
    }
}
