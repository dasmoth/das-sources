package utils;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Collects {
	private Collects() {}
	
	public static String join(Iterable<?> l, String joiner) {
		StringBuilder sb = new StringBuilder();
		Iterator<?> i = l.iterator();
		while (i.hasNext()) {
			sb.append(i.next().toString());
			if (i.hasNext()) {
				sb.append(joiner);
			}
		}
		return sb.toString();
	}
	
	public static String join(Iterable<?> l) {
		return join(l, ",");
	}
	
	public static <K,V> void pushOntoMap(Map<K,List<V>> m, K k, V v)
	{
		List<V> l = m.get(k);
		if (l == null) {
			l = new ArrayList<V>();
			m.put(k, l);
		}
		l.add(v);
	}
	
	public static <K,V> void pushNewOntoMap(Map<K,List<V>> m, K k, V v)
	{
		List<V> l = m.get(k);
		if (l == null) {
			l = new ArrayList<V>();
			m.put(k, l);
		}
		if (!l.contains(v)) {
			l.add(v);
		}
	}
	
	public static <I,O> List<O> map(final Function<I,O> mapper, final List<I> list) {
		return new AbstractList<O>() {
			@Override
			public O get(int i) {
				return mapper.apply(list.get(i));
			}

			@Override
			public int size() {
				return list.size();
			}
		};
	}
	
	public static <I1,I2,O> List<O> map(final Function2<I1,I2,O> mapper, final List<I1> list1, final List<I2> list2) {
		if (list1.size() != list2.size()) {
			throw new IllegalArgumentException("List lengths don't match");
		}
		return new AbstractList<O>() {
			@Override
			public O get(int i) {
				return mapper.apply(list1.get(i), list2.get(i));
			}

			@Override
			public int size() {
				return list1.size();
			}
		};
	}
	
	public static <I1,I2,I3,O> List<O> map(final Function3<I1,I2,I3,O> mapper, final List<I1> list1, final List<I2> list2, final List<I3> list3) {
		if (list1.size() != list2.size() || list1.size() != list3.size()) {
			throw new IllegalArgumentException("List lengths don't match");
		}
		return new AbstractList<O>() {
			@Override
			public O get(int i) {
				return mapper.apply(list1.get(i), list2.get(i), list3.get(i));
			}

			@Override
			public int size() {
				return list1.size();
			}
		};
	}
	
	public static <I> List<I> retainIf(final Function<I,Boolean> test, final List<I> list) {
		List<I> output = new ArrayList<I>();
		for (I i : list) {
			if (test.apply(i)) {
				output.add(i);
			}
		}
		return output;
	}
	
	public static <I, C extends Collection<I>> C  fill(C coll, Iterator<? extends I> iterator)
	{
		while (iterator.hasNext()) {
			coll.add(iterator.next());
		}
		return coll;
	}
	
	public static <I> I reduce(Iterable<? extends I> input, Function2<I, I, I> f)
	{
		Iterator<? extends I> i = input.iterator();
		I cum = i.next();     // Throws NoSuchElementException if input is empty
		while (i.hasNext()) {
			cum = f.apply(cum, i.next());
		}
		return cum;
	}

	public static <X> Set<X> intersection(Set<X> s0, Set<X> s1) {
		if (s0.size() < s1.size()) {
			Set<X> tmp = s1;
			s1 = s0;
			s0 = tmp;
		}
		Set<X> intersect = new HashSet<X>(s1);
		intersect.retainAll(s0);
		return intersect;
	}
	
	/* Far too clever for it's own good.  This is generally going to be less efficient *sigh*
	
	public static <I> List<I> retainIf(final Function<I,Boolean> test, final List<I> list) {
		return new AbstractSequentialList<I>() {

			@Override
			public ListIterator<I> listIterator(int index) {
				return new ListIterator<I>() {

					int publicCursor = 0;
					int cursor = 0;
					int prev = -1;
					int next = -1;
					final int limit = list.size();
					
					
					public void add(I arg0) {
						throw new NotImplementedException();
					}

					public boolean hasNext() {
						if (next < 0) {
							int cc = cursor + 1;
							while (cc < limit && !test.apply(list.get(cc))) {
								++cc;
							}
							if (cc < limit) {
								next = cc;
							}
						}
						return (next >= 0);
					}

					public boolean hasPrevious() {
						if (prev < 0) {
							int cc = cursor - 1;
							while (cc >= 0 && !test.apply(list.get(cc))) {
								--cc;
							}
							if (cc >= 0) {
								next = cc;
							}
						}
						return (next >= 0);
					}

					public I next() {
						if (!hasNext()) {
							throw new NoSuchElementException();
						}
						++publicCursor;
						cursor = next;
						next = -1;
						return list.get(cursor);
					}

					public int nextIndex() {
						return publicCursor + 1;
					}

					public I previous() {
						if (!hasPrevious()) {
							throw new NoSuchElementException();
						}
						--publicCursor;
						cursor = prev;
						prev = -1;
						return list.get(cursor);
					}

					public int previousIndex() {
						return publicCursor - 1; // check!
					}

					public void remove() {
						throw new NotImplementedException();
					}

					public void set(I arg0) {
						throw new NotImplementedException();
					}
					
				};
			}

			@Override
			public int size() {
				int cnt = 0;
				for (I i : this) {
					++cnt;
				}
				return cnt;
			}
			
		};
	}
	
	*/
	
	public static double mean(Collection<? extends Number> l) {
		double tot = 0;
		for (Number n : l) {
			tot += n.doubleValue();
		}
		return tot/l.size();
	}
}
