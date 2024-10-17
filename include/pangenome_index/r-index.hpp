//
// Created by seeskand on 9/21/24.
//
/*
    This is an r-index implementation for grlBWT. These functionalities are needed to be able to do fast locate on the
    BWT created by grlBWT https://github.com/ddiazdom/grlBWT.git
    This implementation is based on the fast_locate implementation in https://github.com/jltsiren/gbwt.git
 */

#ifndef PANGENOME_INDEX_R_INDEX_HPP
#define PANGENOME_INDEX_R_INDEX_HPP

#include <../../deps/grlBWT/include/bwt_io.h>
#include "../../deps/grlBWT/scripts/fm_index.h"
#include <gbwt/algorithms.h>
#include <gbwt/internal.h>

namespace panindexer {

    using namespace gbwt;

/*
  r-index.hpp: A support structure implementing the r-index locate().
*/

//------------------------------------------------------------------------------

/*
  An optional locate() structure based on the r-index. The implementation is based on the
  simplified version in:

    Gagie, Navarro, and Prezza: Fully Functional Suffix Trees and Optimal Text
    Searching in BWT-Runs Bounded Space. Journal of the ACM, 2020.

  We store samples at the beginning of each run instead of at the end, and derive
  SA[i+1] from SA[i]. Because grlBWT is a multi-string BWT, each endmarker is
  considered a separate run.

  The implementation assumes that the product of the number of sequences and the length
  of the longest sequence fits into 64 bits.
*/
    class FastLocate {
    public:
        typedef gbwt::size_type size_type; // Needed for SDSL serialization.

        constexpr static size_type
        NO_POSITION = std::numeric_limits<size_type>::max();

//------------------------------------------------------------------------------

        FastLocate();

        FastLocate(const FastLocate &source);

        FastLocate(FastLocate &&source);

        ~FastLocate();


        explicit FastLocate(const std::string source);

        void swap(FastLocate &another);

        FastLocate &operator=(const FastLocate &source);

        FastLocate &operator=(FastLocate &&source);

        size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const;

        void load(std::istream &in);

        void set_buff_reader(bwt_buff_reader &source) {this->buff_reader = &source;}

        const static std::string EXTENSION; // .ri

//------------------------------------------------------------------------------

        struct Header {
            std::uint32_t tag;
            std::uint32_t version;
            std::uint64_t max_length; // Length of the longest sequence.
            std::uint64_t flags;

            constexpr static std::uint32_t
            TAG = 0x6B3741D8;
            constexpr static std::uint32_t
            VERSION = gbwt::Version::R_INDEX_VERSION;

            Header();

            size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const;

            void load(std::istream &in);

            // Throws `sdsl::simple_sds::InvalidData` if the header is invalid.
            void check() const;

            void setVersion() { this->version = VERSION; }

            void set(std::uint64_t flag) { this->flags |= flag; }

            void unset(std::uint64_t flag) { this->flags &= ~flag; }

            bool get(std::uint64_t flag) const { return (this->flags & flag); }
        };

//------------------------------------------------------------------------------


        std::string grlbwt_file;

        bwt_buff_reader *buff_reader;




        // store the mapping of the symbols to the 0-255 range
        sdsl::int_vector<8> sym_map;

        // store the number of each character in the whole BWT. For example the C[sym_map['A']] gives the total number of As in the whole text
        sdsl::int_vector<64> C;

        // store the size of the whole text
        size_t sequence_size;

        Header header;

        // (sequence id, sequence offset) samples at the start of each run.
        sdsl::int_vector<0> samples;

        // Mark the text positions at the end of a run.
        sdsl::sd_vector<> last;

        // If last[i] = 1, last_to_run[last_rank(i)] is the identifier of the run.
        sdsl::int_vector<0> last_to_run;

        // Run identifier of the first run in each node.
        sdsl::int_vector<0> comp_to_run;

//------------------------------------------------------------------------------

        /*
          Low-level interface: Statistics.
//        */
//
        size_type size() const { return this->samples.size(); }
//        bool empty() const { return (this->size() == 0); }
//
////------------------------------------------------------------------------------
//
//        /*
//          High-level interface. The queries check that the parameters are valid. Iterators
//          must be InputIterators. On error or failed search, the return value will be empty.
//          If the state is non-empty, first will be updated to the packed position corresponding
//          to the first occurrence in the range.
//        */
//
//        SearchState find(node_type node, size_type& first) const;
//
//        template<class Iterator>
//        SearchState find(Iterator begin, Iterator end, size_type& first) const;
//
        gbwt::range_type extend(gbwt::range_type state, size_t sym, size_type& first) const;

//        template<class Iterator>
//        gbwt::range_type extend(gbwt::range_type state, Iterator begin, Iterator end, size_type& first) const;

        std::vector<size_type> locate(gbwt::range_type state, size_type first = NO_POSITION) const;

//        std::vector<size_type> locate(node_type node, gbwt::range_type range, size_type first = NO_POSITION) const
//        {
//            return this->locate(SearchState(node, range), first);
//        }
//
        std::vector<size_type> decompressSA() const;

        std::vector<size_type> decompressDA() const;
//
////------------------------------------------------------------------------------
//
//        /*
//          Low-level interface. The interface assumes that the arguments are valid. This
//          be checked with index->contains(node) and seq_id < index->sequences(). There is
//          no check for the offset.
//        */
//
        size_type pack(size_type seq_id, size_type seq_offset) const
        {
            return seq_id * this->header.max_length + seq_offset;
        }
//
        size_type seqId(size_type offset) const { return offset / this->header.max_length; }
        size_type seqOffset(size_type offset) const { return offset % this->header.max_length; }
//
        std::pair<size_type, size_type> unpack(size_type offset) const
        {
            return std::make_pair(this->seqId(offset), this->seqOffset(offset));
        }



        // This function calculate C and the sym_map from the buff_reader TODO: CHECK THIS
        void calculate_C(){
            this->sequence_size = 0;
            uint8_t buffer[1024]={0};
            size_t sym, freq, k=0;
            size_t sym_freqs[256]={0};
            for(size_t i=0;i<this->buff_reader->size();i++){
                this->buff_reader->read_run(i, sym, freq);
                this->sequence_size += freq;
                for(size_t j=0;j<freq;j++){
                    buffer[k++] = sym;
                    if(k==1024){
                        k=0;
                    }
                }
                sym_freqs[sym]+=freq;
            }

            k=0;
            std::vector<size_t> C_tmp;
            this->sym_map.resize(256);
            sdsl::util::set_to_value(this->sym_map, 0);
            for(size_t i=0;i<256;i++){
                if(sym_freqs[i]!=0){
//                    if(k==0) dummy = i;
                    this->sym_map[i] = k++;
                    C_tmp.push_back(sym_freqs[i]);
                }
            }

            size_t acc=0, tmp;
            for(size_t i=0;i<k;i++){
                tmp = C_tmp[i];
                C_tmp[i] = acc;
                acc+=tmp;
            }
            C_tmp[k] = acc;

            this->C.resize(C_tmp.size());
            for(size_t i=0;i<this->C.size();i++){
                this->C[i] = C_tmp[i];
            }

        }

        size_t tot_strings() const {
            return C[1]-C[0];
        }

        size_type locateFirst() const {
            return this->getSample(0);

        }
//
        size_type locateNext(size_type prev) const;

        gbwt::range_type LF(gbwt::range_type range, size_t sym, bool& starts_with_to, size_t& first_run) const;



        // first is the symbol and the second is the next position index (backtrack)
        std::pair<size_t, size_t> psi(size_t idx);

        size_type total_runs();
//
////------------------------------------------------------------------------------
//
//        /*
//          Internal interface. Do not use.
//        */
//
    private:
        void copy(const FastLocate& source);




        void bwt_index_run_id(unsigned long idx, gbwt::range_type& run, size_type& run_id) ;

        size_t rankAt(size_t pos, size_t symbol, size_t& run_id, size_t& current_position) const;

        size_t rankAt(size_t pos, size_t symbol) const {
            size_t run_id = 0;
            size_t current_position = 0;
            return rankAt(pos, symbol, run_id, current_position);
        };


        size_t bwt_char_at(size_t idx);


//
//        size_type globalRunId(node_type node, size_type run_id) const
//        {
//            return this->comp_to_run[this->index->toComp(node)] + run_id;
//        }
//
        size_type getSample(size_type run_id) const
        {
            return this->samples[run_id];
        }


        size_t get_sequence_size() const { return this->sequence_size; }
//    }; // class FastLocate
//
////------------------------------------------------------------------------------
//
///*
//  Template query implementations.
//*/
//
//    template<class Iterator>
//    SearchState
//    FastLocate::find(Iterator begin, Iterator end, size_type& first) const
//    {
//        if(begin == end) { return SearchState(); }
//
//        SearchState state = this->find(*begin, first);
//        ++begin;
//
//        return this->extend(state, begin, end, first);
//    }
//
//    template<class Iterator>
//    SearchState
//    FastLocate::extend(SearchState state, Iterator begin, Iterator end, size_type& first) const
//    {
//        while(begin != end && !(state.empty()))
//        {
//            state = this->extend(state, *begin, first);
//            ++begin;
//        }
//        return state;
//    }
//
////------------------------------------------------------------------------------
//
    void printStatistics(const FastLocate& index, const std::string& name);
    std::string indexType(const FastLocate&);

//------------------------------------------------------------------------------

    };
}


#endif //PANGENOME_INDEX_R_INDEX_HPP
