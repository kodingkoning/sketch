#include <upcxx/upcxx.hpp>
#include <vector>
#include <set>

using sort_type = int;

size_t message_size_;
size_t domain_size_;
size_t partition_size_;
std::set<sort_type> sorted_set_;
std::vector<std::vector<sort_type>> buckets_;

void sort_handler(std::vector<sort_type> bucket) {
  for (const auto& value : bucket) {
  	sorted_set_.insert(value);
  }
}


void init_sorting(size_t domain_size, size_t message_size = 2048) {
  domain_size_ = domain_size;
  partition_size_ = (domain_size_ + upcxx::rank_n() - 1) / upcxx::rank_n();
  message_size_ = message_size;
  buckets_.resize(upcxx::rank_n());
  upcxx::barrier();
}

void finalize_sorting() {
  for (size_t i = 0; i < buckets_.size(); i++) {
  	auto f = upcxx::rpc(i, sort_handler, buckets_[i]);
  	f.wait();
  }
  upcxx::barrier();
}

void insert_value(sort_type value) {
  size_t dest_proc = value / partition_size_;
  buckets_[dest_proc].push_back(value);

  if (buckets_[dest_proc].size() >= message_size_) {
  	auto f = upcxx::rpc(dest_proc, sort_handler, buckets_[dest_proc]);
  	f.wait();
  	buckets_[dest_proc].clear();
  }
}

int upc_main(int argc, char** argv) {
  upcxx::init();
  size_t num_values = 1000;
  size_t max_value = 10000000;

  // Keep the lowest 100 values.
  size_t num_to_keep = 1500;

  srand48(upcxx::rank_me());

  // Sort values from [0, max_value)
  init_sorting(max_value);

  for (size_t i = 0; i < num_values; i++) {
  	size_t value = lrand48() % max_value;
  	// Sort the value to a remote node.
  	insert_value(value);
  }

  // Send any straggling values.
  finalize_sorting();

  // Get a vector which holds how many values each process has.
  std::vector<size_t> value_count(upcxx::rank_n());
  size_t my_size = sorted_set_.size();
  for (size_t i = 0; i < upcxx::rank_n(); i++) {
  	if (upcxx::rank_me() == i) {
  	  value_count[i] = my_size;
  	}
  	value_count[i] = upcxx::broadcast(value_count[i], i).wait();
  }

  // Get "max_rank" where you can find the num_to_keep
  // smallest values at ranks <= max_rank
  size_t max_rank;
  if (upcxx::rank_me() == 0) {
  	size_t sum = 0;
  	for (size_t i = 0; i < value_count.size(); i++) {
  	  sum += value_count[i];
  	  if (sum >= num_to_keep) {
  	  	max_rank = i;
  	  	break;
  	  }
  	}
  	printf("You can find %lu elements in ranks <= %lu.\n",
  		   num_to_keep, max_rank);
  }

  upcxx::finalize();
  return 0;
}