#ifndef ATLAS_BLAS_FUTURES_H
#define ATLAS_BLAS_FUTURES_H

#include <cassert>
#include <ostream>
extern "C" {
#include <cblas.h>
}

template<typename, typename> class mult_op;
template<typename, typename> class add_op;
template<typename, typename> class sub_op;
template<typename, typename, template<typename, typename> class> class future;
template<typename Left, typename Right> future<Left, Right, mult_op> operator*(Left, Right);
template<typename Left, typename Right> future<Left, Right, add_op> operator+(Left, Right);
template<typename Left, typename Right> future<Left, Right, sub_op> operator-(Left, Right);

template<typename Left, typename Right, template<typename, typename> class Op>
class future {
	public:
		typedef typename Op<Left, Right>::result_type result_type;
		typedef std::size_t size_type;
	private:
		void resolve(result_type&);
		future(Left left, Right right) :
			left(left), right(right) {}
	public:
		Left left;
		Right right;

		size_type cols() { return Op<Left, Right>::cols(left, right); }
		size_type rows() { return Op<Left, Right>::rows(left, right); }
		operator result_type();

	friend std::ostream& operator<<(std::ostream& out, future w)
	{
		out << "(" << w.left << " " << Op<Left, Right>::repr << " " << w.right << ")";
		return out;	
	}

	template<typename, typename> friend class mult_op;
	template<typename, typename> friend class add_op;
	template<typename, typename> friend class sub_op;
};

template<typename Left, typename Right, template<typename, typename> class Op>
void
future<Left, Right, Op>::resolve(result_type& output)
{
	Op<Left, Right>::apply(output, left, right);
}

template<typename Left, typename Right, template<typename, typename> class Op>
future<Left, Right, Op>::operator result_type()
{
	result_type tmp = result_type::make_from(*this);
	resolve(tmp);
	return tmp;
}

template<typename Left, typename Right>
class mult_op {
	public:
		typedef typename mult_op<typename Left::result_type, typename Right::result_type>::result_type result_type;
		typedef std::size_t size_type;
		static const char repr = '*';
	private:
		static void apply(result_type&, Left, Right);
		static future<Left, Right, mult_op> wrap(Left, Right);
		static size_type cols(Left, Right right) { return right.cols(); }
		static size_type rows(Left left, Right) { return left.rows(); }

	template<typename, typename> friend class mult_op;
	friend future<Left, Right, mult_op> operator*<>(Left, Right);
	friend class future<Left, Right, mult_op>;
};

template<typename Left, typename Right>
future<Left, Right, mult_op>
mult_op<Left, Right>::wrap(Left left, Right right)
{
	assert(left.cols() == right.rows());
	return future<Left, Right, mult_op>(left, right);
}

template<typename Left, typename Right>
void
mult_op<Left, Right>::apply(result_type& output, Left left, Right right)
{
	mult_op<typename Left::result_type, typename Right::result_type>::apply(output, left, right);
}

template<typename Right>
class mult_op<double, Right> {
	public:
		typedef typename Right::result_type result_type;
		typedef std::size_t size_type;
		static void apply(result_type&, double, Right);
		static const char repr = '*';
	private:
		static future<double, Right, mult_op> wrap(double, Right);
		static size_type cols(double, Right x) { return x.cols(); }
		static size_type rows(double, Right x) { return x.rows(); }

	template<typename, typename> friend class mult_op;
	friend future<double, Right, mult_op> operator*<>(double, Right);
	friend class future<double, Right, mult_op>;
};

template<typename Right>
future<double, Right, mult_op>
mult_op<double, Right>::wrap(double alpha, Right x)
{
	return future<double, Right, mult_op>(alpha, x);
}

template<typename Right>
void
mult_op<double, Right>::apply(result_type& y, double alpha, Right x)
{
	mult_op<double, typename Right::result_type>::apply(y, alpha, x);
}

template<typename Left, typename Right>
class add_op {
	public:
		typedef typename add_op<typename Left::result_type, typename Right::result_type>::result_type result_type;
		typedef std::size_t size_type;
		static void apply(result_type&, Left, Right); // make private ?
		static const char repr = '+';
	private:
		static future<Left, Right, add_op> wrap(Left, Right);
		static size_type cols(Left, Right right) { return right.cols(); }
		static size_type rows(Left left, Right) { return left.rows(); }

	template<typename, typename> friend class add_op;
	friend future<Left, Right, add_op> operator+<>(Left, Right);
	friend class future<Left, Right, add_op>;
};

template<typename Left, typename Right>
future<Left, Right, add_op>
add_op<Left, Right>::wrap(Left left, Right right)
{
	assert(left.rows() == right.rows());
	return future<Left, Right, add_op>(left, right);
}

template<typename Left, typename Right>
void
add_op<Left, Right>::apply(result_type& output, Left left, Right right)
{
	add_op<typename Left::result_type, typename Right::result_type>::apply(output, left, right);
}

template<typename Left, typename Right>
class sub_op {
	public:
		typedef typename sub_op<typename Left::result_type, typename Right::result_type>::result_type result_type;
		typedef std::size_t size_type;
		static void apply(result_type&, Left, Right); // make private ?
		static const char repr = '-';
	private:
		static future<Left, Right, sub_op> wrap(Left, Right);
		static size_type cols(Left, Right right) { return right.cols(); }
		static size_type rows(Left left, Right) { return left.rows(); }

	template<typename, typename> friend class sub_op;
	friend future<Left, Right, sub_op> operator-<>(Left, Right);
	friend class future<Left, Right, sub_op>;
};

template<typename Left, typename Right>
future<Left, Right, sub_op>
sub_op<Left, Right>::wrap(Left left, Right right)
{
	assert(left.rows() == right.rows());
	return future<Left, Right, sub_op>(left, right);
}

template<typename Left, typename Right>
void
sub_op<Left, Right>::apply(result_type& output, Left left, Right right)
{
	sub_op<typename Left::result_type, typename Right::result_type>::apply(output, left, right);
}

template<typename Left, typename Right>
future<Left, Right, mult_op>
operator*(Left left, Right right)
{
	return mult_op<Left, Right>::wrap(left, right);
}

template<typename Left>
future<double, Left, mult_op>
operator*(Left left, double right)
{
	return right * left;
}

template<typename Left>
future<double, Left, mult_op>
operator/(Left left, double right)
{
	return (1 / right) * left;
}

template<typename Left, typename Right>
future<Left, Right, add_op>
operator+(Left left, Right right)
{
	return add_op<Left, Right>::wrap(left, right);
}

template<typename Left, typename Right>
future<Left, Right, sub_op>
operator-(Left left, Right right)
{
	return sub_op<Left, Right>::wrap(left, right);
}

#endif
