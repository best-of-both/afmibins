#ifndef ATLAS_VECTOR_H
#define ATLAS_VECTOR_H

#include <ostream>

#include "interface.h"
#include "matrix.h"
#include "futures.h"

class dense_matrix;

class vector : public base_matrix {
	public:
		typedef double value_type;
		typedef vector result_type;
		typedef std::size_t index_type;
		typedef std::size_t size_type;
	protected:
		using base_matrix::m_values;
	public:
		value_type& operator[](index_type n) { return m_values[n]; }

		template<typename Wrapper>
		static vector make_from(Wrapper& w) { return vector(w.rows()); }

		template<typename Wrapper>
		vector& operator=(Wrapper);

		vector& operator+=(vector);
		vector& operator-=(vector);

		vector(size_type, double = 0.0);

	/* Unfortunately, we can't do partial specialization here,
	 * but we hope one of the template arguments is vector. */
	template<typename, typename> friend class mult_op;
	template<typename, typename> friend class add_op;
	template<typename, typename> friend class sub_op;

	friend void gemv(vector&, double, dense_matrix, vector, double);
	friend void gemv(vector&, double, sparse_matrix, vector, double);
};

template<>
class mult_op<double, vector> {
	public:
		typedef vector result_type;
		typedef std::size_t size_type;
		static void apply(result_type&, double, vector);
		static const char repr = '*';
	private:
		static future<double, vector, mult_op> wrap(double, vector);
		static size_type cols(double, vector x) { return x.cols(); }
		static size_type rows(double, vector x) { return x.rows(); }

	template<typename, typename> friend class mult_op;
	friend future<double, vector, mult_op> operator*<double, vector>(double, vector);
	friend class future<double, vector, mult_op>;
};

template<>
class add_op<vector, vector> {
	public:
		typedef vector result_type;
		typedef std::size_t size_type;
		static void apply(result_type&, vector, vector); // make private ?
		static const char repr = '+';
	private:
		static future<vector, vector, add_op> wrap(vector, vector);
		static size_type cols(vector, vector right) { return right.cols(); }
		static size_type rows(vector, vector right) { return right.rows(); }

	template<typename, typename> friend class add_op;
	friend future<vector, vector, add_op> operator+<>(vector, vector);
	friend class future<vector, vector, add_op>;
};

template<>
class sub_op<vector, vector> {
	public:
		typedef vector result_type;
		typedef std::size_t size_type;
		static void apply(result_type&, vector, vector); // make private ?
		static const char repr = '-';
	private:
		static future<vector, vector, sub_op> wrap(vector, vector);
		static size_type cols(vector, vector right) { return right.cols(); }
		static size_type rows(vector, vector right) { return right.rows(); }

	template<typename, typename> friend class sub_op;
	friend future<vector, vector, sub_op> operator-<>(vector, vector);
	friend class future<vector, vector, sub_op>;
};

std::ostream& operator<<(std::ostream&, vector&);

#endif
