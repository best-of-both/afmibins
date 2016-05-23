#ifndef ATLAS_SPARSE_MATRIX_H
#define ATLAS_SPARSE_MATRIX_H

#include "interface.h"
#include "futures.h"
#include "vector.h"
#include "matrix.h"

class sparse_matrix : public base_matrix {
	public:
		typedef sparse_matrix result_type;
		using base_matrix::index_type;
		using base_matrix::size_type;
	protected:
		size_type m_length;
		index_type* m_row_starts;
		index_type* m_cols;
		using base_matrix::m_values;

		size_type num_elements() { return m_length; }
	public:
		template<typename Wrapper>
		static sparse_matrix make_from(Wrapper& w) { return sparse_matrix(w.rows(), w.cols()); }

		sparse_matrix(size_type, size_type);
		~sparse_matrix();

	template<typename, typename> friend class mult_op;
	template<typename, typename> friend class add_op;
	template<typename, typename> friend class sub_op;

	friend void gemv(vector&, double, sparse_matrix, vector, double);
};

template<>
class mult_op<sparse_matrix, vector> {
	public:
		typedef vector result_type;
		typedef std::size_t size_type;
		static const char repr = '*';
	private:
		static void apply(result_type&, sparse_matrix, vector); // make private ?
		static future<sparse_matrix, vector, mult_op> wrap(sparse_matrix, vector);
		static size_type cols(sparse_matrix, vector x) { return x.cols(); }
		static size_type rows(sparse_matrix a, vector) { return a.rows(); }

	template<typename, typename> friend class mult_op;
	friend future<sparse_matrix, vector, mult_op> operator*<>(sparse_matrix, vector);
	friend class future<sparse_matrix, vector, mult_op>;
};

template<>
class mult_op<double, sparse_matrix> {
	public:
		typedef sparse_matrix result_type;
		typedef std::size_t size_type;
		static void apply(result_type&, double, sparse_matrix);
		static const char repr = '*';
	private:
		static future<double, sparse_matrix, mult_op> wrap(double, sparse_matrix);
		static size_type cols(double, sparse_matrix x) { return x.cols(); }
		static size_type rows(double, sparse_matrix x) { return x.rows(); }

	template<typename, typename> friend class mult_op;
	friend future<double, sparse_matrix, mult_op> operator*<double, sparse_matrix>(double, sparse_matrix);
	friend class future<double, sparse_matrix, mult_op>;
};

template<>
class add_op<sparse_matrix, sparse_matrix> {
	public:
		typedef sparse_matrix result_type;
		typedef std::size_t size_type;
		static void apply(result_type&, sparse_matrix, sparse_matrix); // make private ?
		static const char repr = '+';
	private:
		static future<sparse_matrix, sparse_matrix, add_op> wrap(sparse_matrix, sparse_matrix);
		static size_type cols(sparse_matrix, sparse_matrix right) { return right.cols(); }
		static size_type rows(sparse_matrix, sparse_matrix right) { return right.rows(); }

	template<typename, typename> friend class add_op;
	friend future<sparse_matrix, sparse_matrix, add_op> operator+<>(sparse_matrix, sparse_matrix);
	friend class future<sparse_matrix, sparse_matrix, add_op>;
};

template<>
class sub_op<sparse_matrix, sparse_matrix> {
	public:
		typedef sparse_matrix result_type;
		typedef std::size_t size_type;
		static void apply(result_type&, sparse_matrix, sparse_matrix); // make private ?
		static const char repr = '-';
	private:
		static future<sparse_matrix, sparse_matrix, sub_op> wrap(sparse_matrix, sparse_matrix);
		static size_type cols(sparse_matrix, sparse_matrix right) { return right.cols(); }
		static size_type rows(sparse_matrix, sparse_matrix right) { return right.rows(); }

	template<typename, typename> friend class sub_op;
	friend future<sparse_matrix, sparse_matrix, sub_op> operator-<>(sparse_matrix, sparse_matrix);
	friend class future<sparse_matrix, sparse_matrix, sub_op>;
};

std::ostream& operator<<(std::ostream&, sparse_matrix&);

#endif
