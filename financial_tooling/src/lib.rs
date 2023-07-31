use plotly::common::{Mode, Title};
use plotly::layout::{Axis, Layout};
use plotly::{Bar, Plot, Scatter};
use pyo3::prelude::*;
use std::f64::consts::E;

#[pyclass]
#[derive(Clone, Debug)]
pub struct PWLinearInterpolator {
    #[pyo3(get, set)]
    pub abcissae: Vec<f64>,
    #[pyo3(get, set)]
    pub ordinates: Vec<f64>,
    #[pyo3(get, set)]
    pub length: usize,
}

#[pymethods]
impl PWLinearInterpolator {
    #[new]
    pub fn new() -> PWLinearInterpolator {
        PWLinearInterpolator {
            abcissae: Vec::new(),
            ordinates: Vec::new(),
            length: 0,
        }
    }
    pub fn from(&mut self, abcissae: Vec<f64>, ordinates: Vec<f64>) -> PWLinearInterpolator {
        let mut temp_vec: Vec<f64> = abcissae.clone();
        let length: usize = abcissae.len();
        temp_vec.sort_by(|a, b| a.partial_cmp(b).unwrap());

        if abcissae.len() != ordinates.len() {
            println!("Length of abcissae and ordinates not the same.",);
            PWLinearInterpolator {
                abcissae,
                ordinates,
                length,
            }
        } else if temp_vec != abcissae {
            println!("Warning: xvalues out of order");
            PWLinearInterpolator {
                abcissae,
                ordinates,
                length,
            }
        } else {
            PWLinearInterpolator {
                abcissae,
                ordinates,
                length,
            }
        }
    }

    pub fn extend(&mut self, xval: f64, yval: f64) {
        self.abcissae.push(xval);
        self.abcissae.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let index: usize = self.abcissae.iter().position(|&r| r == xval).unwrap();
        self.ordinates.insert(index, yval);
        self.length += 1;
    }

    pub fn eval(&mut self, x: f64) -> f64 {
        let mut boolvect: Vec<bool> = Vec::new();
        for xval in self.abcissae.clone() {
            if x <= xval {
                boolvect.push(true);
            } else {
                boolvect.push(false);
            }
        }

        let index: usize = boolvect
            .iter()
            .position(|&r| r == true)
            .unwrap_or(self.length);

        if self.length == 0 {
            println!("Empty Interpolation");
            return 0.0;
        }
        if index == 0 {
            return self.ordinates[0];
        } else {
            if x > *self.abcissae.last().unwrap() {
                // Out of range extrapolation
                println!("Warning: Extrapolating out of range!");
                return self.ordinates[index - 1];
            } else {
                let alpha: f64 = (x - self.abcissae[index - 1])
                    / (self.abcissae[index] - self.abcissae[index - 1]);
                let value: f64 =
                    (1.0 - alpha) * self.ordinates[index - 1] + alpha * self.ordinates[index];
                return value;
            }
        }
    }

    pub fn delta(&mut self, x: f64, bump_index: usize) -> f64 {
        let mut boolvect: Vec<bool> = Vec::new();
        for xval in self.abcissae.clone() {
            if x <= xval {
                boolvect.push(true);
            } else {
                boolvect.push(false);
            }
        }

        let index: usize = boolvect.iter().position(|&r| r == true).unwrap_or(0);

        let bump_index_i32 = bump_index as i32;

        if bump_index_i32 < (index as i32 - 1) || bump_index > index {
            return 0.0;
        }

        if self.length == 0 {
            println!("Empty Interpolation");
            return 0.0; // need to figure out error returning
        } else if self.length == 1 {
            return 1.0;
        } else {
            if x > *self.abcissae.last().unwrap() {
                // Out of range extrapolation
                println!("Warning: Extrapolating out of range!");
                if bump_index == (self.length - 1) {
                    return 1.0;
                } else {
                    return 0.0;
                }
            } else if x == self.abcissae[0] {
                return self.ordinates[0];
            } else {
                let alpha: f64 = (x - self.abcissae[index - 1])
                    / (self.abcissae[index] - self.abcissae[index - 1]);
                if bump_index == (index - 1) {
                    return 1.0 - alpha;
                }
                if bump_index == index {
                    return alpha;
                } else {
                    return 99.0;
                }
            }
        }
    }
}

#[pyfunction]
fn solver(interpolator: PWLinearInterpolator, bond: &AbstractBond) -> f64 {
    let tolerance: f64 = 0.000001;

    let dates = bond.dates.clone();
    let coupons = bond.coupons.clone();
    let price = bond.market_price.clone();
    let maturity = bond.maturity.clone();

    let f = |eta: f64| {
        let mut _interpolator = interpolator.clone();
        _interpolator.extend(maturity, eta);

        let mut sum = -price;

        for i in 0..dates.len() {
            sum += E.powf(-1.0 * _interpolator.eval(dates[i]) * dates[i]) * coupons[i];
        }
        return sum;
    };
    let fprime = |eta: f64| {
        let mut _interpolator = interpolator.clone();
        _interpolator.extend(maturity, eta);

        let max_index = _interpolator.length - 1;
        let mut sum = 0.0;

        for i in 0..dates.len() {
            println!("{:?}", _interpolator.delta(dates[i], max_index));
            sum += _interpolator.delta(dates[i], max_index)
                * dates[i]
                * E.powf(-1.0 * _interpolator.eval(dates[i]) * dates[i])
                * coupons[i];
        }
        return -sum;
    };

    let mut old_eta = 0.0;
    let mut new_eta = old_eta - (f(old_eta) / fprime(old_eta));

    while (new_eta - old_eta).abs() > tolerance {
        old_eta = new_eta;
        new_eta = old_eta - (f(old_eta) / fprime(old_eta));
    }

    return new_eta;
}

#[pyclass]
#[derive(FromPyObject, Debug)]
struct AbstractBond {
    #[pyo3(get, set)]
    pub coupons: Vec<f64>,
    #[pyo3(get, set)]
    pub dates: Vec<f64>,
    #[pyo3(get, set)]
    pub maturity: f64,
    #[pyo3(get, set)]
    pub length: usize,
    #[pyo3(get, set)]
    pub market_price: f64,
}

#[pymethods]
impl AbstractBond {
    #[new]
    pub fn new() -> AbstractBond {
        AbstractBond {
            coupons: Vec::new(),
            dates: Vec::new(),
            length: 0,
            maturity: 0.0,
            market_price: 0.0,
        }
    }
    pub fn build(&mut self, rates: Vec<f64>, times: Vec<f64>, price: f64) -> AbstractBond {
        let l: usize = rates.len();
        let maturity = times.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
        if rates.len() == 0 {
            println!("No coupons given.")
        } else if rates.len() != times.len() {
            println!("Coupon and Tensor vectors are not of equal length.")
        }
        AbstractBond {
            coupons: rates,
            dates: times.clone(),
            maturity: *maturity,
            length: l,
            market_price: price,
        }
    }
    pub fn build_from_coupon_bond(
        &mut self,
        maturity: f64,
        face: f64,
        rate: f64,
        frequency: f64,
    ) -> AbstractBond {
        let mut bond = AbstractBond::new();
        if frequency == 0.0 {
            let paydates = vec![maturity];
            let payments = vec![face];
            return bond.build(payments, paydates, 0.0);
        } else {
            let mut paydates: Vec<f64> = Vec::new();
            let mut payments: Vec<f64> = Vec::new();
            let period = 1.0 / frequency;
            let mut date = maturity;
            while date > 0.0 {
                paydates.push(date);
                date -= period;
            }

            paydates.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let num_of_coupons = paydates.len();
            let coupon = (rate / 100.0) * face / frequency;

            for i in 0..num_of_coupons {
                payments.push(coupon);
            }

            payments[num_of_coupons - 1] += face;

            return bond.build(payments, paydates, 0.0);
        }
    }
    pub fn plot_payments(&mut self) {
        let layout = Layout::new()
            .title(Title::new("Coupon Bond Payments"))
            .y_axis(Axis::new().title(Title::new("Payments")))
            .x_axis(Axis::new().title(Title::new("Payment Dates (years)")));
        let mut plot = Plot::new();
        let t = Bar::new(self.dates.clone(), self.coupons.clone());
        plot.add_trace(t);
        plot.set_layout(layout);
        plot.show();
    }
    pub fn price(&mut self, mut curve: Curve, date: f64) -> f64 {
        let mut pv: f64 = 0.0;

        for i in 0..self.length {
            if self.dates[i] < date {
                continue;
            }
            pv += curve.discount_factor(self.dates[i] - date) * self.coupons[i];
        }
        return pv;
    }
    pub fn ytm(&mut self, price: f64) -> f64 {
        let mut price = price;
        if price == 0.0 {
            price = self.market_price;
        }
        if price == 0.0 {
            println!("No price available!")
        }
        let tolerance = 0.000001;

        let frequency = (1.0 / (self.dates[1] - self.dates[0])).floor();

        let f = |y: f64| {
            let mut sum =
                self.coupons[0] / (1.0 + (y / frequency)).powf(frequency * self.dates[0] + 1.0);

            for j in 1..self.length {
                sum +=
                    self.coupons[j] / (1.0 + (y / frequency)).powf(frequency * self.dates[j] + 1.0)
            }

            return sum - price;
        };
        let fprime = |y: f64| {
            let mut sum = self.dates[0] * self.coupons[0]
                / (1.0 + (y / frequency)).powf(frequency * self.dates[0] + 1.0);

            for j in 1..self.length {
                sum += self.dates[j] * self.coupons[j]
                    / (1.0 + (y / frequency)).powf(frequency * self.dates[j] + 1.0)
            }

            return -sum;
        };

        let mut old_y = 0.0;
        let mut new_y = old_y - f(old_y) / fprime(old_y);

        while (new_y - old_y).abs() > tolerance {
            old_y = new_y;
            new_y = old_y - f(old_y) / fprime(old_y);
        }

        return 100.0 * new_y;
    }
}

#[pyclass]
#[derive(Clone)]
struct Curve {
    #[pyo3(get, set)]
    pub tenors: Vec<f64>,
    #[pyo3(get, set)]
    pub rates: Vec<f64>,
    #[pyo3(get, set)]
    pub length: usize,
    #[pyo3(get, set)]
    pub interpolator: Option<PWLinearInterpolator>,
}

#[pymethods]
impl Curve {
    #[new]
    pub fn new() -> Curve {
        Curve {
            tenors: Vec::new(),
            rates: Vec::new(),
            length: 0,
            interpolator: None,
        }
    }
    pub fn build_from_rates(&mut self, dates: Vec<f64>, rates: Vec<f64>) -> Curve {
        if dates.len() != rates.len() {
            println!("Rate and tenor vectors are of unequal length.")
        }
        let length = dates.len();
        let mut i = PWLinearInterpolator::new();
        let interpolator = Some(i.from(dates.clone(), rates.clone()));

        Curve {
            tenors: dates,
            rates,
            length,
            interpolator,
        }
    }
    pub fn build_from_bonds(&mut self, bond_list: Vec<AbstractBond>) {
        let mut bond_list = bond_list;

        bond_list.sort_by(|a, b| a.maturity.partial_cmp(&b.maturity).unwrap());
        self.interpolator = Some(PWLinearInterpolator::new());

        for bond in bond_list.iter() {
            let new_rate = solver(self.interpolator.clone().unwrap(), bond);
            self.interpolator
                .as_mut()
                .unwrap()
                .extend(bond.maturity, new_rate);
            self.tenors.push(bond.maturity);
            self.rates.push(new_rate);
            self.length += 1;
        }
    }
    pub fn shift(&mut self, delta: f64) -> Curve {
        let delta = delta / 100.0;
        let mut new_rates = self.rates.clone();
        let new_tenors = self.tenors.clone();

        for i in 0..self.tenors.len() {
            new_rates[i] = self.rates[i] + delta;
        }

        return self.build_from_rates(new_tenors, new_rates);
    }
    pub fn bump(&mut self, delta: f64, index: usize) -> Curve {
        let delta = delta / 100.0;
        let mut new_rates = self.tenors.clone();
        let new_tenors = self.tenors.clone();
        new_rates[index] = new_rates[index] + delta;
        return self.build_from_rates(new_tenors, new_rates);
    }
    pub fn get_yield(&mut self, time: f64) -> f64 {
        return 100.0 * self.interpolator.as_mut().unwrap().eval(time);
    }
    pub fn discount_factor(&mut self, time: f64) -> f64 {
        return E.powf(self.interpolator.as_mut().unwrap().eval(time) * time);
    }
    pub fn spot_rate(&mut self, time: f64, compounding: usize) -> f64 {
        if compounding == 0 {
            return 100.0 * (1.0 / self.discount_factor(time) - 1.0) / time;
        } else {
            let df = self.discount_factor(time);
            return 100.0 * compounding as f64 * (df.powf(-1.0 / compounding as f64 / time) - 1.0);
        }
    }
    pub fn plot_yields(&mut self, max_tenor: f64) {
        let mut max_tenor = max_tenor;
        if max_tenor <= 0.0 {
            max_tenor = *self.tenors.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
        }

        let mut t: Vec<f64> = Vec::new();
        for i in (0..(max_tenor * 100.0) as usize).step_by(10) {
            t.push(i as f64 * 0.01);
        }

        let mut y: Vec<f64> = Vec::new();
        for tenor in t.iter() {
            let yld = self.get_yield(*tenor);
            y.push(yld);
        }

        let layout = Layout::new()
            .title(Title::new("Yields"))
            .y_axis(Axis::new().title(Title::new("yield")))
            .x_axis(Axis::new().title(Title::new("tenor")));
        let mut plot = Plot::new();
        let t = Scatter::new(t.clone(), y).mode(Mode::LinesMarkers);
        plot.add_trace(t);
        plot.set_layout(layout);
        plot.show();
    }
    pub fn plot_discount_factors(&mut self, max_tenor: f64) {
        let mut max_tenor = max_tenor;
        if max_tenor <= 0.0 {
            max_tenor = *self.tenors.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
        }

        let mut t: Vec<f64> = Vec::new();
        for i in (0..(max_tenor * 10.0) as usize + 1).step_by(1) {
            t.push(i as f64 / 10.0);
        }

        let mut y: Vec<f64> = Vec::new();
        for tenor in t.iter() {
            let yld = self.discount_factor(*tenor);
            y.push(yld);
        }

        let layout = Layout::new()
            .title(Title::new("Discount Factors"))
            .y_axis(Axis::new().title(Title::new("discount factor")))
            .x_axis(Axis::new().title(Title::new("tenor")));
        let mut plot = Plot::new();
        let t = Scatter::new(t.clone(), y).mode(Mode::LinesMarkers);
        plot.add_trace(t);
        plot.set_layout(layout);
        plot.show();
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn financial_tooling(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PWLinearInterpolator>()?;
    m.add_class::<AbstractBond>()?;
    m.add_class::<Curve>()?;
    m.add_function(wrap_pyfunction!(solver, m)?)?;
    Ok(())
}
