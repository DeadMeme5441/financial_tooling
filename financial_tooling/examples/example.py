import financial_tooling

# a = [1.0, 2.0, 3.0, 4.0, 5.0]
# b = [10.0, 12.0, 14.0, 16.0, 18.0]
# interp = financial_tooling.PWLinearInterpolator(a, b)

# print(interp.eval(3.0))

bond = financial_tooling.AbstractBond()

bond1 = bond.build_from_coupon_bond(1.0, 10000.0, 7.0, 2)
bond1.market_price = 10050.0
bond2 = bond.build_from_coupon_bond(1.5, 50000.0, 4.0, 2)
bond2.market_price = 46500.0
bond3 = bond.build_from_coupon_bond(2.0, 100000.0, 5.0, 2)
bond3.market_price = 96500.0
bond4 = bond.build_from_coupon_bond(2.5, 100000.0, 7.0, 2)
bond4.market_price = 98000.0
# bond1.plot_payments()

bonds = [bond1, bond2, bond3, bond4]

yc = financial_tooling.Curve()
yc.build_from_bonds(bonds)
print(yc.rates)
yc.plot_yields(0.0)
