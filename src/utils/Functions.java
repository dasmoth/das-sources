package utils;

public class Functions {
	private Functions() {
	}
	
	// Ugh, why can't the compiler do some sensible type-inference here?
	
	public static <X> Function<X,X> identity(Class<X> dummyClass) {
		return new Function<X, X>() {
			public X apply(X param1) {
				return param1;
			}
		};
	}
	
	public static <I,O,M> Function<I,O> compose(final Function<I,M> inner, final Function<M,O> outer) {
		return new Function<I, O>() {
			public O apply(I param1) {
				return outer.apply(inner.apply(param1));
			}
		};
	}
}
